#include "LidDrivenCavity.h"
#include "cblas.h"
#include <iostream>
#include <iomanip>


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

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
    delete[] omegaInterior;
    delete[] omegaRight;
    delete[] omegaLeft;
    delete[] omegaTop;
    delete[] omegaBottom;
    
    delete[] psiInterior;
    delete[] psiRight;
    delete[] psiLeft;
    delete[] psiTop;
    delete[] psiBottom;
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
}

void LidDrivenCavity::Initialise()
{
    // Calculating resuable constant and setting it
    udy = 2*U/dy;
    
    // Initialise as zeros array for initial conditions
    // t = 0, omega(x,y) = 0, psi(x,y) = 0
    omegaInterior = new double[interiorNarr]{};
    omegaRight = new double[interiorNy]{};
    omegaLeft = new double[interiorNy]{};
    omegaTop = new double[interiorNx]{};
    omegaBottom = new double[interiorNx]{};
    
    psiInterior = new double[interiorNarr]{};
    psiRight = new double[interiorNy]{};
    psiLeft = new double[interiorNy]{};
    psiTop = new double[interiorNx]{};
    psiBottom = new double[interiorNx]{};
    
    // Generating laplacian matrix
    SetSymmetricBandedLaplacianMatrix();
}

void LidDrivenCavity::SetSymmetricBandedLaplacianMatrix()
{
    // Initialising zeros matrix
    symmetricBandedLaplacianMatrix = new double[interiorNarr*(interiorNy+1)]();
    
    double alpha = 1 / dx / dx;
    double beta = 1 / dy / dy;
    double gamma = 2*(alpha+beta);
    
    for (unsigned int i=0; i<interiorNarr; i++) {
        // Column-major storage of laplacian banded matrix
        symmetricBandedLaplacianMatrix[i*(interiorNy+1) + 0] = gamma;
        if ((i+1)%interiorNy!=0) {
            symmetricBandedLaplacianMatrix[i*(interiorNy+1) + 1] = -1*beta;
        }
        if (i<interiorNy*(interiorNx-1)) {
            symmetricBandedLaplacianMatrix[i*(interiorNy+1) + interiorNy] = -1*alpha;
        }
    }
//    PrintMatrix(symmetricBandedLaplacianMatrix,interiorNy+1,interiorNarr,false);
}


void LidDrivenCavity::SetVorticityBoundaryConditions()
{
    // Setting top bottom boundary conditions
    for (unsigned int i=0; i<interiorNx; i++) {
        // Setting the top BC
        omegaTop[i] = (psiTop[i]-psiInterior[i+(interiorNy-2)*interiorNx])*2/dy/dy-udy;
        // Setting the bottom BC
        omegaBottom[i] = (psiBottom[i]-psiInterior[i+(1)*interiorNx])*2/dy/dy;
    }
    
    // Setting left right boundary conditions
    for (unsigned int i=0; i<Ny; i++) {
        // Setting the left BC
        omegaLeft[i] = (psiLeft[i]-psiInterior[i*interiorNx+1])*2/dx/dx;
        // Setting the right BC
        omegaRight[i] = (psiRight[i]-psiInterior[i*interiorNx+(interiorNx-2)])*2/dx/dx;
    }
}

void LidDrivenCavity::SetInteriorVorticity()
{
    for (unsigned int i=0; i<interiorNx; i++) {
        for (unsigned int j=0; j<interiorNy; j++) {
            psiInterior[i*interiorNy + j] = rand()%100;
        }
        psiRight[i] = rand()%100;
        psiLeft[i] = rand()%100;
        psiTop[i] = rand()%100;
        psiBottom[i] = rand()%100;
    }
    PrintMatrix(psiInterior,interiorNy,interiorNx,false);
    PrintMatrix(psiRight,interiorNy,1,false);
    PrintMatrix(psiLeft,interiorNy,1,false);
    PrintMatrix(psiTop,1,interiorNx,false);
    PrintMatrix(psiBottom,1,interiorNx,false);

    cblas_dsbmv(CblasColMajor, CblasLower, interiorNarr, interiorNy, 1.0, symmetricBandedLaplacianMatrix, interiorNy+1, psiInterior, 1, 0.0, omegaInterior, 1);
    cblas_daxpy(interiorNy, -1/dx/dx, psiLeft, 1, omegaInterior, 1); // Left
    cblas_daxpy(interiorNy, -1/dy/dy, psiBottom, 1, omegaInterior, interiorNy); // Bottom
    cblas_daxpy(interiorNy, -1/dx/dx, psiRight, 1, omegaInterior+((interiorNx-1)*interiorNy), 1); // Right
    cblas_daxpy(interiorNy, -1/dy/dy, psiTop, 1, omegaInterior+(interiorNy-1), interiorNy); // Top
}

void LidDrivenCavity::UpdateInteriorVorticity()
{
    double* term3 = new double[interiorNarr]();
    
    cblas_dsbmv(CblasColMajor, CblasLower, interiorNarr, interiorNy, 1.0, symmetricBandedLaplacianMatrix, interiorNy+1, omegaInterior, 1, 0.0, term3, 1);
    cblas_daxpy(interiorNy, -1/dx/dx, omegaLeft, 1, term3, 1); // Left
    cblas_daxpy(interiorNy, -1/dy/dy, omegaBottom, 1, term3, interiorNy); // Bottom
    cblas_daxpy(interiorNy, -1/dx/dx, omegaRight, 1, term3+((interiorNx-1)*interiorNy), 1); // Right
    cblas_daxpy(interiorNy, -1/dy/dy, omegaTop, 1, term3+(interiorNy-1), interiorNy); // Top
    
    cblas_daxpy(interiorNarr, -dt/Re, term3, 1, omegaInterior, 1);
    PrintMatrix(omegaInterior,interiorNy,interiorNx,false);
    delete[] term3;
}

void LidDrivenCavity::SolvePoissonProblem()
{
    
}

void LidDrivenCavity::Integrate()
{
    SetVorticityBoundaryConditions();
    SetInteriorVorticity();
    UpdateInteriorVorticity();
    SolvePoissonProblem();
}

// REMOVE(?)
void LidDrivenCavity::BruteForceSetInteriorVorticity()
{
//    for (unsigned int i=1; i<(Nx-1); i++) {
//        for (unsigned int j=1; j<(Ny-1); j++) {
//            omega[i*Nx+j] = -(psi[i*Nx+(j+1)]-2*psi[i*Nx+(j)]+psi[i*Nx+(j-1)]) /dx/dx - (psi[(i+1)*Nx+j]-2*psi[i*Nx+j]+psi[(i-1)*Nx+j]) /dy/dy;
//        }
//    }
}

void LidDrivenCavity::BruteForceUpdateInteriorVorticity()
{
//    for (unsigned int i=1; i<(Nx-1); i++) {
//        for (unsigned int j=1; j<(Ny-1); j++) {
//            double term1 = (psi[i*Nx+(j+1)]-psi[i*Nx+(j-1)])*(omega[(i+1)*Nx+j]-omega[(i-1)*Nx+j]);
//            double term2 = (omega[i*Nx+(j+1)]-omega[i*Nx+(j-1)])*(psi[(i+1)*Nx+j]-psi[(i-1)*Nx+j]);
//            double term3 = (omega[i*Nx+(j+1)]-2*omega[i*Nx+(j)]+omega[i*Nx+(j-1)]) /dx/dx;
//            double term4 = (omega[(i+1)*Nx+j]-2*omega[i*Nx+j]+omega[(i-1)*Nx+j]) /dy/dy;
//
//            omega[i*Nx+j] = dt/4/dx/dy * (term1 - term2) + dt/Re * (term3 + term4);
//        }
//    }
//    PrintOmegaMatrix();
}

void LidDrivenCavity::BruteForceIntegrate()
{
    BruteForceSetInteriorVorticity();
    BruteForceUpdateInteriorVorticity();
    SolvePoissonProblem();
}



// Remove these
void LidDrivenCavity::PrintOmegaMatrix() {
//    for (int x=0; x<Nx; x++) {
//        for (int y=0; y<Ny; y++) {
//            cout << setprecision(3) <<setw(7) << left << omega[x*Nx+y];
//        }
//        cout << endl;
//    }
}
