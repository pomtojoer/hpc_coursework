/**
    High Performane Computing Coursework
    LidDrivenCavity.cpp
    Purpose: Defines the class member functions for the Poisson2DSolver class. In general, it implements a serial and parallel solver for a 2D Poisson Equation (Dirichlet Problem) using Lapack and Scalapack.

    @author Sean Chai
    @version 1.0 23/03/20
*/

#include "LidDrivenCavity.h"
#include "Poisson2DSolver.h"
#include "cblas.h"
#include "mpi.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>


// Defines the constructor for the LidDrivenCavity class
LidDrivenCavity::LidDrivenCavity()
{
}


// Defines the destructor for the LidDrivenCavity class
LidDrivenCavity::~LidDrivenCavity()
{
    // Deallocating memory of class arrays
    delete[] w;
    delete[] s;
}


/**
    Configure and initialises MPI as well as segmenting the w and s array into evenly distributed chunks for each process.
    They are segmented in the order of their storage position.
    
    @return void
*/
void LidDrivenCavity::SetMPIConfig()
{
    // Checking if MPI was initially initalised
    int MPIInitialised;
    MPI_Initialized(&MPIInitialised);
    if (!MPIInitialised){
        cout << "Error: MPI was not initialised." << endl;
        throw std::exception();
    } else {
        // Get the rank and comm size on each process.
        MPI_Comm_rank(MPI_COMM_WORLD, &MPIRank);
        MPI_Comm_size(MPI_COMM_WORLD, &MPISize);
        
        // Checking if serial or parallel
        if (MPISize > 1) {
            // Calculating the number of elemets to allocate to each bin (partition)
            int spacing = (int) floor((double) interiorNarr / (double)MPISize);
            int remainder = (int)interiorNarr%MPISize;
            
            int tag = 0;
            int source = 0;
            
            // Generating the coordinate pairs and allocating to the bins from the root
            if (MPIRank == source) {
                // Generating i-j inner coordinate pairs
                int* iCoord = new int[interiorNarr];
                int* jCoord = new int[interiorNarr];
                for (unsigned int i=0; i<interiorNx; i++) {
                    for (unsigned int j=0; j<interiorNy; j++) {
                        iCoord[i*interiorNy+j] = i+1;
                        jCoord[i*interiorNy+j] = j+1;
                    }
                }
                
                // Calculating the start and end points of the bins
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
                
                // Sending the array of binned i,j coordinates to the different processes
                for (int dest=0; dest<MPISize; dest++) {
                    int startPos = startingPositions[dest];
                    int endPos = endingPositions[dest];
                    int lenArr = endPos-startPos+1;
                    // splitting the inner coordinates of the inner matrix into respective bins
                    if (dest==0) {
                        coordArrLen = lenArr;
                        iInnerCoords = new int[coordArrLen];
                        jInnerCoords = new int[coordArrLen];
                        std::copy(iCoord+startPos, iCoord+endPos+1, iInnerCoords);
                        std::copy(jCoord+startPos, jCoord+endPos+1, jInnerCoords);
                    } else {
                        // Creating the array containing the local chunk of coordinates
                        int* tempiCoords = new int[lenArr];
                        int* tempjCoords = new int[lenArr];
                        std::copy(iCoord+startPos, iCoord+endPos+1, tempiCoords);
                        std::copy(jCoord+startPos, jCoord+endPos+1, tempjCoords);
                        
                        //Sending the array containing the coordinates to the individual processes
                        MPI_Send(&lenArr, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
                        MPI_Send(tempiCoords, lenArr, MPI_INT, dest, tag, MPI_COMM_WORLD);
                        MPI_Send(tempjCoords, lenArr, MPI_INT, dest, tag, MPI_COMM_WORLD);
                        delete[] tempiCoords;
                        delete[] tempjCoords;
                    }
                }
                
                // Deallocatin the temporary arrays
                delete[] startingPositions;
                delete[] endingPositions;
                
                delete[] iCoord;
                delete[] jCoord;
            } else {
                // Receiving the binned i,j coordinates from the root
                MPI_Recv(&coordArrLen, 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                iInnerCoords = new int[coordArrLen];
                jInnerCoords = new int[coordArrLen];
                MPI_Recv(iInnerCoords, coordArrLen, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(jInnerCoords, coordArrLen, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
}


/**
    Function to set the length of the domain in xy directions.
    
    @param xlen The length of the domain in the x direction
    @param ylen The length of the domain in the y direction
    @return void
*/
void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    Lx = xlen;
    Ly = ylen;
}


/**
    Function to set the grid size in the xy domain. The function also calculates and sets other useful variables.
    
    @param nx The number of grid points in the x direction
    @param ny The number of grid points in the y direction
    @return void
*/
void LidDrivenCavity::SetGridSize(unsigned int nx, unsigned int ny)
{
    Nx = nx;
    Ny = ny;
    
    // Total matrix size
    narr = Nx*Ny;
    
    // Inner matrix dimensions
    interiorNx = Nx-2;
    interiorNy = Ny-2;
    
    // Inner matrix size
    interiorNarr = interiorNx * interiorNy;
}


/**
    Function to set the partition size in the xy domain.
 
    @param nx The number of grid points in the x direction
    @param ny The number of grid points in the y direction
    @return void
*/
void LidDrivenCavity::SetPartitionSize(unsigned int px, unsigned int py)
{
    Px = px;
    Py = py;
}


/**
    Function to set the time step
 
    @param deltat The time step of each iteration
    @return void
*/
void LidDrivenCavity::SetTimeStep(double deltat)
{
    dt = deltat;
}


/**
    Function to set the final time
 
    @param finalt The final time
    @return void
*/
void LidDrivenCavity::SetFinalTime(double finalt)
{
    T = finalt;
}


/**
    Function to set the reynolds number
 
    @param re The reynolds number
    @return void
*/
void LidDrivenCavity::SetReynoldsNumber(double re)
{
    Re = re;
}


/**
    Function to set the grid spacing
 
    @param deltax The grid spacing in the x direction
    @param deltay The grid spacing in the y direction
    @return void
*/
void LidDrivenCavity::SetGridSpacing(double deltax, double deltay)
{
    // Function to set the grid spacing
    dx = deltax;
    dy = deltay;
}


/**
    Function to initialise the w and s arrays
 
    @param
    @return void
*/
void LidDrivenCavity::Initialise()
{
    // Initialise as zeros array for initial conditions
    // t = 0, omega(x,y) = 0, psi(x,y) = 0
    w = new double[narr]{};
    s = new double[narr]{};
}


/**
    Function to set the vorticity matrix (w) boundary conditions
 
    @param
    @return void
*/
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


/**
    Function to set the interior vorticity (w) at time t. Implements a second order central difference scheme for the spacial derivatives
 
    @param
    @return void
*/
void LidDrivenCavity::SetInteriorVorticity()
{
    if (MPISize > 1) {
        // Parallel version
        // Initialising matrix to store local portion of the matrix
        double* tempW = new double[narr]{};
        
        // Calculating the interior vorticity for the assigned segment on the given process
        for (int k=0; k<coordArrLen; k++) {
            int i = iInnerCoords[k];
            int j = jInnerCoords[k];
            
            tempW[i*Ny+j] = -(s[i*Ny+(j+1)]-2*s[i*Ny+(j)]+s[i*Ny+(j-1)]) /dy/dy - (s[(i+1)*Ny+j]-2*s[i*Ny+j]+s[(i-1)*Ny+j]) /dx/dx;
        }

        if (MPIRank>0) {
            // sending the calculated local w matrix to the root process
            MPI_Send(tempW, narr, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        } else {
            // receiving each individual local segment and combining them on root to the global
            for (int src=1; src<MPISize; src++) {
                double* tempWseg = new double[narr]{};
                MPI_Recv(tempWseg, narr, MPI_DOUBLE, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                cblas_daxpy(narr, 1.0, tempWseg, 1, tempW, 1);
                delete[] tempWseg;
            }
        }
        
        // updating the w matrix
        delete[] w;
        w = tempW;
    } else {
        // Serial version
        // Calculating the interior vorticity for whole matrix w
        for (unsigned int i=1; i<(Nx-1); i++) {
            for (unsigned int j=1; j<(Ny-1); j++) {
                w[i*Ny+j] = -(s[i*Ny+(j+1)]-2*s[i*Ny+(j)]+s[i*Ny+(j-1)]) /dy/dy - (s[(i+1)*Ny+j]-2*s[i*Ny+j]+s[(i-1)*Ny+j]) /dx/dx;
            }
        }
    }
    
    
}


/**
    Function to set the interior vorticity (w) at time t+dt. Implements a second order central difference scheme for the spacial derivatives.
 
    @param
    @return void
*/
void LidDrivenCavity::UpdateInteriorVorticity()
{
    // Initialising matrix to store updated w. W is not updated immediately as the neighbouring segments
    // rely on the t values to calculate t+dt
    double* temp = new double[narr]{};
    
    if (MPISize > 1) {
        // Parallel version
        // Broadcasting the original w matrix at time t to all processes
        MPI_Bcast(w, narr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        // Calculating the updated interior vorticity for the assigned segment on the given process
        for (int k=0; k<coordArrLen; k++) {
            int i = iInnerCoords[k];
            int j = jInnerCoords[k];
            double term1 = (s[(i+1)*Ny+j]-s[(i-1)*Ny+j])*(w[i*Ny+(j+1)]-w[i*Ny+(j-1)]);
            double term2 = (s[i*Ny+(j+1)]-s[i*Ny+(j-1)])*(w[(i+1)*Ny+j]-w[(i-1)*Ny+j]);
            double term3 = (w[i*Ny+(j+1)]-2*w[i*Ny+(j)]+w[i*Ny+(j-1)]) /dy/dy;
            double term4 = (w[(i+1)*Ny+j]-2*w[i*Ny+j]+w[(i-1)*Ny+j]) /dx/dx;
            
            temp[i*Ny+j] = dt/4/dx/dy * (term1 - term2) + dt/Re * (term3 + term4);
        }
        
        if (MPIRank>0) {
            // sending the calculated local w matrix to the root process
            MPI_Send(temp, narr, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        } else {
            // receiving each individual local segment and combining them on root to the global
            for (int src=1; src<MPISize; src++) {
                double* tempW = new double[narr];
                MPI_Recv(tempW, narr, MPI_DOUBLE, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                cblas_daxpy(narr, 1.0, tempW, 1, temp, 1);
                delete[] tempW;
            }
        }
        
    } else {
        // Serial version
        // Calculating the updated interior vorticity for whole matrix w
        for (unsigned int i=1; i<(Nx-1); i++) {
            for (unsigned int j=1; j<(Ny-1); j++) {
                double term1 = (s[(i+1)*Ny+j]-s[(i-1)*Ny+j])*(w[i*Ny+(j+1)]-w[i*Ny+(j-1)]);
                double term2 = (s[i*Ny+(j+1)]-s[i*Ny+(j-1)])*(w[(i+1)*Ny+j]-w[(i-1)*Ny+j]);
                double term3 = (w[i*Ny+(j+1)]-2*w[i*Ny+(j)]+w[i*Ny+(j-1)]) /dy/dy;
                double term4 = (w[(i+1)*Ny+j]-2*w[i*Ny+j]+w[(i-1)*Ny+j]) /dx/dx;
                
                temp[i*Ny+j] = dt/4/dx/dy * (term1 - term2) + dt/Re * (term3 + term4);
            }
        }
    }
    
    // updating the current w (at time t) matrix with the new w matrix (at time t+dt)
    cblas_daxpy(narr,1,temp,1,w,1);
    
    // deallocating temporary matrix
    delete[] temp;
}


/**
    Function that integrates all the steps together to solve for the streamfunction and vorticity iteratively
 
    @param
    @return void
*/
void LidDrivenCavity::Integrate()
{
    // calculating the maindiagonal and superdiagonal values of the A matrix
    double alpha = 1 / dx / dx;
    double beta = 1 / dy / dy;
    double gamma = 2*(alpha+beta);
    
    // Setting up the poisson solver
    Poisson2DSolver* poissonSolver = new Poisson2DSolver();
    poissonSolver->SetVariables((int)Nx,(int)Ny, alpha, beta, gamma);
    
    if (MPISize > 1) {
        // Parallel version
        if (MPIRank==0) {
            // Generating the Laplacian matrix on the root process
            poissonSolver->GenerateScalapackMatrixAHat();
            scalapackMatrix = poissonSolver->GetScalapackMatrixAHat();
            scalapackMatrixNx = poissonSolver->GetScalapackMatrixAHatNx();
            scalapackMatrixNy = poissonSolver->GetScalapackMatrixAHatNy();
        }
        // Sending the laplacian matrix to the other processes (root)
        // Receiving and setting the laplacian matrix from the root process (other processes)
        MPI_Bcast(&scalapackMatrixNx, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&scalapackMatrixNy, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (MPIRank > 0) {
            scalapackMatrix = new double[scalapackMatrixNx*scalapackMatrixNy];
        }
        MPI_Bcast(scalapackMatrix, scalapackMatrixNx * scalapackMatrixNy, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (MPIRank > 0) {
            poissonSolver->SetScalapackMatrixAHat(scalapackMatrix, scalapackMatrixNx, scalapackMatrixNy);
        }
        
        // Initialising BLACS for ScaLAPACK
        poissonSolver->InitialiseScalapack(Px,Py);
        // Prefactoring the A matrix for faster solve
        poissonSolver->PrefactorMatrixAHatParallel();
    } else {
        // Serial version
        // Generating the Laplacian matrix
        poissonSolver->GenerateLapackMatrixAHat();
        // Prefactoring the A matrix for faster solve
        poissonSolver->PrefactorMatrixAHatSerial();
    }
    
    // Initialising a variable to keep track of the current time during the solve
    double tnow = 0.0;
    do {
        // Output to the terminal the current percentage completion
        if (MPIRank==0) {
            string strout = "\r" + to_string((int)ceil((tnow/T*100))) + "% completed";
            cout << strout;
        }
        
//        if (MPIRank==0) {
//            cout << "At initial" << endl;
//            cout << "omega" << endl;
//            PrintMatrix(w,Ny,Nx,false);
//            cout << "psi" << endl;
//            PrintMatrix(s,Ny,Nx,false);
//        }

        SetInteriorVorticity(); // Setting interior vorticity
//        if (MPIRank==0) {
//            cout << "After setting vorticity at t" << endl;
//            cout << "omega" << endl;
//            PrintMatrix(w,Ny,Nx,false);
//            cout << "psi" << endl;
//            PrintMatrix(s,Ny,Nx,false);
//        }

        SetVorticityBoundaryConditions();   // Setting vorticity boundary conditions
//        if (MPIRank==0) {
//            cout << "After setting BC" << endl;
//            cout << "omega" << endl;
//            PrintMatrix(w,Ny,Nx,false);
//            cout << "psi" << endl;
//            PrintMatrix(s,Ny,Nx,false);
//        }
        
        UpdateInteriorVorticity();  // Updating interior vorticity
//        if (MPIRank==0) {
//            cout << "After setting vorticity at t+dt" << endl;
//            cout << "omega" << endl;
//            PrintMatrix(w,Ny,Nx,false);
//            cout << "psi" << endl;
//            PrintMatrix(s,Ny,Nx,false);
//        }
        
        
        // Broadcasting the assembled vorticity matrix to all other processes
        if (MPISize > 1) {
            MPI_Bcast(w, narr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        
        // Setting the streamfunction and vorticity matrices for the poisson solver to use
        poissonSolver->SetVectors(s, w);
        
        // Solving based on serial or parallel
        if (MPISize > 1) {
            poissonSolver->SolveParallel();
        } else {
            poissonSolver->SolveSerial();
        }
        
        // Updating the current streamfunction (at time t) with the new stream function (at time t+dt)
        if (MPIRank==0) {
            poissonSolver->Updatex(s);
        }
        
        // Broadcasting the updated streamfunction to all other processes
        if (MPISize > 1) {
            MPI_Bcast(s, narr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        
//        if (MPIRank==0) {
//            cout << "After solve" << endl;
//            cout << "omega" << endl;
//            PrintMatrix(w,Ny,Nx,false);
//            cout << "psi" << endl;
//            PrintMatrix(s,Ny,Nx,false);
//        }
        
        // Updating the current time
        tnow += dt;
    } while (tnow < T);
}


/**
    Function that outputs the streamfunction and vorticity matrices into a .txt file
 
    @param
    @return void
*/
void LidDrivenCavity::GeneratePlotData()
{
    if (MPIRank==0) {
        // Creating filename from parameters used
        string filename = "Data/" +  to_string((int)Lx) + "_" + to_string((int)Ly) + "_" + to_string(Nx) + "_" + to_string(Ny) + "_" + to_string((int)Re) + "_data.txt";
        
        // Opening file to write. If present before, overwrite it.
        ofstream outputfile(filename, ios::out | ios::trunc);
        if (outputfile.is_open()) {
            // Output the initial parameters
            outputfile << "Lx, Ly, Nx, Ny, dt, T, Re" << endl;
            outputfile << Lx << "," << Ly <<  "," << Nx <<  "," << Ny <<  "," << dt <<  "," << T <<  "," << Re << endl;
            
            // Output of the vorticity
            outputfile << endl << "wdata" << endl;
            for (unsigned int j=0; j<Ny; j++) {
                for (unsigned int i=0; i<Nx; i++) {
                    outputfile << w[j+i*Ny] << "    ";
                }
                outputfile << endl;
            }
            
            // Output of the streamfunction
            outputfile << endl << "sdata" << endl;
            for (unsigned int j=0; j<Ny; j++) {
                for (unsigned int i=0; i<Nx; i++) {
                    outputfile << s[j+i*Ny] << "    ";
                }
                outputfile << endl;
            }
        }
    }
}


// Helper functions
//void PrintMatrix(double* mat, int rows, int cols, bool isRowMajor) {
//    for (int i=0; i<rows; i++) {
//        for (int j=0; j<cols; j++) {
//            if (isRowMajor == false) {
//                cout << setprecision(5) << setw(10) << left << mat[i+j*rows];
//            } else {
//                cout << setprecision(5) << setw(10) << left << mat[i*cols+j];
//            }
//        }
//        cout << endl;
//    }
//    cout << endl;
//}
//
//void printVector(double* vec, int n) {
//    cout << "The vector looks like: " << endl;
//    for (int i=0; i<n; i++) {
//        cout << vec[i] << endl;
//    }
//}
//
//void printVector(int* vec, int n) {
//    cout << "The vector looks like: " << endl;
//    for (int i=0; i<n; i++) {
//        cout << vec[i] << endl;
//    }
//    cout << endl;
//}

