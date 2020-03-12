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

void printVector(double* vec, int n) {
    cout << "The vector looks like: " << endl;
    for (int i=0; i<n; i++) {
        cout << vec[i] << endl;
    }
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
    
    psiInterior = new double[interiorNarr]{7,49,73,78,23,9,87,3,27};
    psiRight = new double[interiorNy]{58,40,29};
    psiLeft = new double[interiorNy]{30,65,40};
    psiTop = new double[interiorNx]{72,92,12};
    psiBottom = new double[interiorNx]{44,42,3};
    
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
}


void LidDrivenCavity::SetVorticityBoundaryConditions()
{
    // Setting top bottom boundary conditions
    for (unsigned int i=0; i<interiorNx; i++) {
        // Setting the top BC
        omegaTop[i] = (psiTop[i]-psiInterior[(interiorNy-1)+i*interiorNy])*2/dy/dy-udy;
        // Setting the bottom BC
        omegaBottom[i] = (psiBottom[i]-psiInterior[(0)+i*interiorNy])*2/dy/dy;
    }
    
    // Setting left right boundary conditions
    for (unsigned int i=0; i<interiorNy; i++) {
        // Setting the left BC
        omegaLeft[i] = (psiLeft[i]-psiInterior[i])*2/dx/dx;
        // Setting the right BC
        omegaRight[i] = (psiRight[i]-psiInterior[i+(interiorNx-1)*interiorNy])*2/dx/dx;
    }
}

void LidDrivenCavity::SetInteriorVorticity()
{
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
    double* t1 = new double[interiorNarr]{};
    double* t4 = new double[interiorNarr]{};
    double* term1 = new double[interiorNarr]{};
    double* term3 = new double[interiorNarr]{};
    
    // implementation of t1
    // (psi(i+1,j) - psi(i-1,j))
    cblas_dcopy(interiorNarr-interiorNy, psiInterior+interiorNy, 1, t1, 1);
    cblas_daxpy(interiorNarr-interiorNy, -1.0, psiInterior, 1, t1+interiorNy, 1);
    cblas_daxpy(interiorNy, -1.0, psiLeft, 1, t1, 1);
    cblas_daxpy(interiorNy, 1.0, psiRight, 1, t1+(interiorNx-1)*interiorNy, 1);
    
    // implementation of t4
    // (omega(i+1,j) - omega(i-1,j))
    cblas_dcopy(interiorNarr-interiorNy, omegaInterior+interiorNy, 1, t4, 1);
    cblas_daxpy(interiorNarr-interiorNy, -1.0, omegaInterior, 1, t4+interiorNy, 1);
    cblas_daxpy(interiorNy, -1.0, omegaLeft, 1, t4, 1);
    cblas_daxpy(interiorNy, 1.0, omegaRight, 1, t4+(interiorNx-1)*interiorNy, 1);

    // implementation of term 1 (LHS of interior vorticity at t+dt)
    for (int i=0; i<interiorNx; i++) {
        for (int j=0; j<interiorNy; j++) {
            if (j==0) {
                term1[j+i*interiorNy] = (t1[j+i*interiorNy]*(omegaInterior[(j+1)+i*interiorNx]-omegaBottom[i])) - (t4[j+i*interiorNy]*(psiInterior[(j+1)+i*interiorNx]-psiBottom[i]));
            } else if (j==interiorNy-1) {
                term1[j+i*interiorNy] = (t1[j+i*interiorNy]*(omegaTop[i]-omegaInterior[(j-1)+i*interiorNx])) - (t4[j+i*interiorNy]*(psiTop[i]-psiInterior[(j-1)+i*interiorNx]));
            } else {
                term1[j+i*interiorNy] = (t1[j+i*interiorNy]*(omegaInterior[(j+1)+i*interiorNx]-omegaInterior[(j-1)+i*interiorNx])) - (t4[j+i*interiorNy]*(psiInterior[(j+1)+i*interiorNx]-psiInterior[(j-1)+i*interiorNx]));
            }
        }
    }
    
    // Implementing RHS of interior vorticity at t+dt
    cblas_dsbmv(CblasColMajor, CblasLower, interiorNarr, interiorNy, 1.0, symmetricBandedLaplacianMatrix, interiorNy+1, omegaInterior, 1, 0.0, term3, 1);
    cblas_daxpy(interiorNy, -1/dx/dx, omegaLeft, 1, term3, 1); // Left
    cblas_daxpy(interiorNy, -1/dy/dy, omegaBottom, 1, term3, interiorNy); // Bottom
    cblas_daxpy(interiorNy, -1/dx/dx, omegaRight, 1, term3+((interiorNx-1)*interiorNy), 1); // Right
    cblas_daxpy(interiorNy, -1/dy/dy, omegaTop, 1, term3+(interiorNy-1), interiorNy); // Top
    
    // Updating omega with new values
    cblas_daxpy(interiorNarr, -dt/Re, term3, 1, omegaInterior, 1);
    cblas_daxpy(interiorNarr, dt/4/dx/dy, term1, 1, omegaInterior, 1);
    
    // Deallocating memory of temp variables
    delete[] t1;
    delete[] t4;
    delete[] term1;
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



// Remove these
void LidDrivenCavity::PrintOmegaMatrix() {
//    for (int x=0; x<Nx; x++) {
//        for (int y=0; y<Ny; y++) {
//            cout << setprecision(3) <<setw(7) << left << omega[x*Nx+y];
//        }
//        cout << endl;
//    }
}
