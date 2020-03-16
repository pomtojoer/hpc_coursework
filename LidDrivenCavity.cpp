#include "LidDrivenCavity.h"
#include "Poisson2DSolver.h"
#include "cblas.h"
#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <cmath>

void PrintMatrix(double* mat, int rows, int cols, bool isRowMajor) {
    for (int i=0; i<rows; i++) {
        for (int j=0; j<cols; j++) {
            if (isRowMajor == false) {
                cout << setw(10) << left << mat[i+j*rows];
            } else {
                cout << setw(10) << left << mat[i*cols+j];
            }
        }
        cout << endl;
    }
    cout << endl;
}

void printVector(double* vec, int n) {
    cout << "The vector looks like: " << endl;
    for (int i=0; i<n; i++) {
        cout << vec[i] << endl;
    }
}

void printVector(int* vec, int n) {
    cout << "The vector looks like: " << endl;
    for (int i=0; i<n; i++) {
        cout << vec[i] << endl;
    }
    cout << endl;
}

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
    delete[] w;
    delete[] s;
}

void LidDrivenCavity::SetMPIConfig()
{
    MPI_Initialized(&MPIInitialised);
    if (!MPIInitialised){
        cout << "Error: MPI was not initialisde." << endl;
        throw std::exception();
    } else {
        // Get the rank and comm size on each process.
        MPI_Comm_rank(MPI_COMM_WORLD, &MPIRank);
        MPI_Comm_size(MPI_COMM_WORLD, &MPISize);
        
        if (MPISize > 1) {
            int spacing = (int) floor((double) interiorNarr / (double)MPISize);
            int remainder = (int)interiorNarr%MPISize;
            
            int tag = 0;
            int source = 0;
            
            if (MPIRank == source) {
                // Generating i-j inner coordinate pairs
                
                int* iCoord = new int[interiorNarr];
                int* jCoord = new int[interiorNarr];
                for (unsigned int i=0; i<interiorNx; i++) {
                    for (unsigned int j=0; j<interiorNx; j++) {
                        iCoord[i*interiorNy+j] = i+1;
                        jCoord[i*interiorNy+j] = j+1;
                    }
                }
                
                int* startingPositions = new int [MPISize];
                int* endingPositions = new int [MPISize];
                
                int currentLocation = 0;
                for (int i=0; i<MPISize; i++) {
                    startingPositions[i] = currentLocation;
                    currentLocation += spacing;
                    if (remainder>0) {
                        currentLocation++;
                        remainder--;
                    }
                    endingPositions[i] = currentLocation-1;
                }
                
                for (int dest=0; dest<MPISize; dest++) {
                    int startPos = startingPositions[dest];
                    int endPos = endingPositions[dest];
                    int lenArr = endPos-startPos+1;
                    if (dest==0) {
                        coordArrLen = lenArr;
                        iInnerCoords = new int[coordArrLen];
                        jInnerCoords = new int[coordArrLen];
                        std::copy(iCoord+startPos, iCoord+endPos+1, iInnerCoords);
                        std::copy(jCoord+startPos, jCoord+endPos+1, jInnerCoords);
                    } else {
                        int* tempiCoords = new int[lenArr];
                        int* tempjCoords = new int[lenArr];
                        std::copy(iCoord+startPos, iCoord+endPos+1, tempiCoords);
                        std::copy(jCoord+startPos, jCoord+endPos+1, tempjCoords);
                        MPI_Send(&lenArr, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
                        MPI_Send(tempiCoords, lenArr, MPI_INT, dest, tag, MPI_COMM_WORLD);
                        MPI_Send(tempjCoords, lenArr, MPI_INT, dest, tag, MPI_COMM_WORLD);
                        delete[] tempiCoords;
                        delete[] tempjCoords;
                    }
                }
                delete[] startingPositions;
                delete[] endingPositions;
                
                delete[] iCoord;
                delete[] jCoord;
            }
            else {
               MPI_Recv(&coordArrLen, 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               iInnerCoords = new int[coordArrLen];
               jInnerCoords = new int[coordArrLen];
               MPI_Recv(iInnerCoords, coordArrLen, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(jInnerCoords, coordArrLen, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    Lx = xlen;
    Ly = ylen;
}

void LidDrivenCavity::SetGridSize(unsigned int nx, unsigned int ny)
{
    Nx = nx;
    Ny = ny;
    
    narr = Nx*Ny;
    
    interiorNx = Nx-2;
    interiorNy = Ny-2;
    
    interiorNarr = interiorNx * interiorNy;
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
    dt = deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
    T = finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
    Re = re;
}

void LidDrivenCavity::SetGridSpacing(double deltax, double deltay)
{
    dx = deltax;
    dy = deltay;
    
    alpha = 1 / dx / dx;
    beta = 1 / dy / dy;
    gamma = 2*(alpha+beta);
}

void LidDrivenCavity::Initialise()
{
    // Initialise as zeros array for initial conditions
    // t = 0, omega(x,y) = 0, psi(x,y) = 0
    w = new double[narr]{};
    s = new double[narr]{7,49,73,58,30,72,44,78,23,9,40,65,92,42,87,3,27,29,40,12,3,69,9,57,60};
}


void LidDrivenCavity::SetVorticityBoundaryConditions()
{
    if(MPIRank==0) {
        // Setting top bottom boundary conditions
        for (unsigned int i=0; i<Nx; i++) {
            // Setting the top BC
            w[i*Ny+(Ny-1)] = (s[i*Ny+(Ny-1)]-s[i*Ny+(Ny-2)]) * 2/dy/dy - 2*U/dy;
            // Setting the bottom BC
            w[i*Ny] = (s[i*Ny]-s[i*Ny+1]) * 2/dy/dy;
        }

        // Setting left right boundary conditions
        for (unsigned int i=0; i<Ny; i++) {
            // Setting the left BC
            w[i] = (s[i]-s[i+Ny]) * 2/dx/dx;
            // Setting the right BC
            w[i+Ny*(Nx-1)] = (s[i+Ny*(Nx-1)]-s[i+Ny*(Nx-2)]) * 2/dx/dx;
        }
    }
}

void LidDrivenCavity::SetInteriorVorticity()
{
    if (MPISize > 1) {
//        cout << "Rank " << MPIRank << " Start: " << startingPosition << " End: " << endingPosition << endl;
        for (int k=0; k<coordArrLen; k++) {
            int i = iInnerCoords[k];
            int j = jInnerCoords[k];
            
            w[i*Nx+j] = -(s[i*Nx+(j+1)]-2*s[i*Nx+(j)]+s[i*Nx+(j-1)]) /dx/dx - (s[(i+1)*Nx+j]-2*s[i*Nx+j]+s[(i-1)*Nx+j]) /dy/dy;
        }
        
        if (MPIRank>0) {
            MPI_Send(w, narr, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        } else {
            for (int src=1; src<MPISize; src++) {
                double* tempW = new double[narr];
                MPI_Recv(tempW, narr, MPI_DOUBLE, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                cblas_daxpy(narr, 1.0, tempW, 1, w, 1);
                delete[] tempW;
            }
        }
    } else {
        for (unsigned int i=1; i<(Nx-1); i++) {
            for (unsigned int j=1; j<(Ny-1); j++) {
                w[i*Ny+j] = -(s[i*Ny+(j+1)]-2*s[i*Ny+(j)]+s[i*Ny+(j-1)]) /dx/dx - (s[(i+1)*Ny+j]-2*s[i*Ny+j]+s[(i-1)*Ny+j]) /dy/dy;
            }
        }
    }
    
    
}

void LidDrivenCavity::UpdateInteriorVorticity()
{
    double* temp = new double[narr]{};
    
    for (unsigned int i=1; i<(Nx-1); i++) {
        for (unsigned int j=1; j<(Ny-1); j++) {
            double term1 = (s[(i+1)*Ny+j]-s[(i-1)*Ny+j])*(w[i*Ny+(j+1)]-w[i*Ny+(j-1)]);
            double term2 = (s[i*Ny+(j+1)]-s[i*Ny+(j-1)])*(w[(i+1)*Ny+j]-w[(i-1)*Ny+j]);
            double term3 = (w[i*Ny+(j+1)]-2*w[i*Ny+(j)]+w[i*Ny+(j-1)]) /dx/dx;
            double term4 = (w[(i+1)*Ny+j]-2*w[i*Ny+j]+w[(i-1)*Ny+j]) /dy/dy;
            
            temp[i*Nx+j] = dt/4/dx/dy * (term1 - term2) + dt/Re * (term3 + term4);
        }
    }
    
    cblas_daxpy(narr,1,temp,1,w,1);
    delete[] temp;
}

void LidDrivenCavity::Integrate()
{
    double* news = new double[narr];
    
    // Note that in our case we do not need to impose the boundary conditions on the w matrix
    // due to the fact that the edges of the s matrix are always zeros (stream function zero)
//    Poisson2DSolver* poissonSolver = new Poisson2DSolver();
//    poissonSolver->Initialise(news, w, Nx, Ny);
//    poissonSolver->InitialiseMPI();
//    poissonSolver->GenerateScalapackMatrixAHat(alpha,beta,gamma);
    
    if (MPIRank==0) {
        cout << "At initial" << endl;
        cout << "omega" << endl;
        PrintMatrix(w,Ny,Nx,false);
        cout << "psi" << endl;
        PrintMatrix(s,Ny,Nx,false);
    }

    SetVorticityBoundaryConditions();
    if (MPIRank==0) {
        cout << "After setting BC" << endl;
        cout << "omega" << endl;
        PrintMatrix(w,Ny,Nx,false);
        cout << "psi" << endl;
        PrintMatrix(s,Ny,Nx,false);
    }

    SetInteriorVorticity();
    if (MPIRank==0) {
        cout << "After setting vorticity at t" << endl;
        cout << "omega" << endl;
        PrintMatrix(w,Ny,Nx,false);
        cout << "psi" << endl;
        PrintMatrix(s,Ny,Nx,false);
    }
//    cout << "omega" << endl;
//    PrintMatrix(w,Ny,Nx,false);
//    cout << "psi" << endl;
//    PrintMatrix(s,Ny,Nx,false);
//
//    UpdateInteriorVorticity();
//    cout << "After setting vorticity at t+dt" << endl;
//    cout << "omega" << endl;
//    PrintMatrix(w,Ny,Nx,false);
//    cout << "psi" << endl;
//    PrintMatrix(s,Ny,Nx,false);
    
//    poissonSolver->ApplyBoundaryConditions();
//
//    cout << "Testing started here" << endl;
//    printVector(psiTop,interiorNx);
//    printVector(psiBottom,interiorNx);
//    printVector(psiLeft,interiorNy);
//    printVector(psiRight,interiorNy);
//    PrintMatrix(psiInterior,interiorNy,interiorNx,false);
//    cout << endl;
//    printVector(omegaTop,interiorNx);
//    printVector(omegaBottom,interiorNx);
//    printVector(omegaLeft,interiorNy);
//    printVector(omegaRight,interiorNy);
//    PrintMatrix(omegaInterior,interiorNy,interiorNx,false);
    
    if (MPISize>1) {
        if (MPIRank>0) {
            int source = 0;
            MPI_Recv(s, narr, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            for (int dest=1; dest<MPISize; dest++) {
                MPI_Send(s, narr, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            }
        }
    }
}

void LidDrivenCavity::GeneratePlots()
{
    
}
