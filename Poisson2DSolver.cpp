#include "Poisson2DSolver.h"
#include "cblas.h"
#include <iostream>
#include <iomanip>

Poisson2DSolver::Poisson2DSolver()
{
}

Poisson2DSolver::~Poisson2DSolver()
{
    delete[] AHat;
    delete[] xHat;
    delete[] b;
    
    delete[] topBC;
    delete[] bottomBC;
    delete[] leftBC;
    delete[] rightBC;
}

void PrintMatrix2(double* mat, int rows, int cols, bool isRowMajor) {
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

void printVector2(double* vec, int n) {
    cout << "The vector looks like: " << endl;
    for (int i=0; i<n; i++) {
        cout << vec[i] << endl;
    }
}

void Poisson2DSolver::SetVariables(int nx, int ny) {
    Nx = nx;
    Ny = ny;
    
    bHatNx = Nx-2;
    bHatNy = Ny-2;
}

void Poisson2DSolver::SetVectors(double* xhat, double* bVec)
{
    xHat = xhat;
    b = bVec;
    
    cblas_dcopy(bHatNy,b+1,1,leftBC,1);             // Left BC
    cblas_dcopy(bHatNy,b+(Nx-1)*Ny+1,1,rightBC,1);  // Right BC
    cblas_dcopy(bHatNy,b+(1)*Ny+1,Ny,topBC,1);      // Top
    cblas_dcopy(bHatNy,b+(2)*Ny-1,Ny,bottomBC,1);   // Bottom BC
    
    printVector2(leftBC,bHatNy);
    printVector2(rightBC,bHatNy);
    printVector2(topBC,bHatNy);
    printVector2(bottomBC,bHatNy);
}

void Poisson2DSolver::InitialiseMPI()
{
}

void Poisson2DSolver::GenerateScalapackMatrixAHat(double alpha, double beta, double gamma)
{
    aHatNx = bHatNx*bHatNy;
    aHatNy = (bHatNy*2+1);
    AHat = new double[aHatNx*aHatNy]();
    for (int i=0; i<bHatNy*bHatNx; i++) {
        /* Column-major storage of laplacian banded symmetric matrix for interiorNx=3, interiorNy=3
         | *  *  *  *  *  *  *  *  * |
         | *  *  *  *  *  *  *  *  * |
         | *  *  *  *  *  *  *  *  * |
         | g  g  g  g  g  g  g  g  g |
         |-b -b  * -b -b  * -b -b  * |
         | *  *  *  *  *  *  *  *  * |
         |-a -a -a -a -a -a  *  *  * |
         
         Where the repeating cell is of size interiorNy*(interiorNy+1) and is repeated interiorNx times
         */
        
        AHat[i*(bHatNy*2+1) + bHatNy] = gamma;
        if ((i+1)%bHatNy!=0) {
            AHat[i*(bHatNy*2+1) + 1 + bHatNy] = -1*beta;
        }
        if (i<bHatNy*(bHatNx-1)) {
            AHat[i*(bHatNy*2+1) + 2*bHatNy] = -1*alpha;
        }
    }
    PrintMatrix2(AHat,aHatNy,aHatNx,false);
}

double* Poisson2DSolver::GetScalapackMatrixAHat() {
    return AHat;
}

int Poisson2DSolver::GetScalapackMatrixAHatNx() {
    return aHatNx;
}

int Poisson2DSolver::GetScalapackMatrixAHatNy() {
    return aHatNy;
}

void Poisson2DSolver::SetScalapackMatrixAHat(double* ahat, int ahatnx, int ahatny) {
    AHat = ahat;
    aHatNx = ahatnx;
    aHatNy = ahatny;
}
