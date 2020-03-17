#include "Poisson2DSolver.h"
#include "cblas.h"
#include "mpi.h"
#include <iostream>
#include <iomanip>

#define F77NAME(x) x##_
extern "C" {
    // For serial solve
    void F77NAME(dgbtrf) (const int& m, const int& n, const int& kl,
                          const int& ku, double* A, const int& lda,
                          int* ipiv, int& info);

    void F77NAME(dgbtrs) (const char& trans, const int& n, const int& kl,
                          const int &ku, const int& nrhs, const double* A,
                          const int& lda, const int* ipiv, double* b,
                          const int& ldb, int& info);
    

    void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                        const int& lda, int * ipiv, double * B,
                        const int& ldb, int& info);


}

Poisson2DSolver::Poisson2DSolver()
{
}

Poisson2DSolver::~Poisson2DSolver()
{
    delete[] AHat;
    delete[] xHat;
    delete[] bHat;
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

void Poisson2DSolver::SetVariables(int nx, int ny, double alphaVar, double betaVar, double gammaVar) {
    Nx = nx;
    Ny = ny;
    
    alpha = alphaVar;
    beta = betaVar;
    gamma =gammaVar;
    
    bHatNx = Nx-2;
    bHatNy = Ny-2;
    
    aHatNx = bHatNx*bHatNy;
    aHatNy = (bHatNy*4+1);
    
    bHat = new double[bHatNx*bHatNy];
    xHat = new double[bHatNx*bHatNy];
    AHat = new double[aHatNx*aHatNy];
}

void Poisson2DSolver::SetVectors(double* xVec, double* bVec)
{
    for (int i=0; i<bHatNx; i++){
        for (int j=0; j<bHatNy; j++) {
            bHat[i*bHatNy+j] = bVec[(i+1)*Ny+(j+1)];
        }
    }
    
    cblas_daxpy(bHatNy, alpha, xVec+1, 1, bHat, 1);                              // Left BC
    cblas_daxpy(bHatNy, alpha, xVec+(Nx-1)*Ny+1, 1, bHat+bHatNy*(bHatNx-1), 1);  // Right BC
    cblas_daxpy(bHatNx, beta, xVec+(2)*Ny-1, Ny, bHat+(bHatNy-1), bHatNy);       // Top BC
    cblas_daxpy(bHatNx, beta, xVec+(1)*Ny, Ny, bHat, bHatNy);                    // Bottom BC
}

void Poisson2DSolver::InitialiseMPI()
{
}

void Poisson2DSolver::GenerateScalapackMatrixAHat()
{
    for (int i=0; i<aHatNx; i++) {
        /* Column-major storage of laplacian banded symmetric matrix for interiorNx=3, interiorNy=3
         | *  *  *  *  *  *  *  *  * |
         | *  *  *  *  *  *  *  *  * |
         | *  *  *  *  *  *  *  *  * |
         | *  *  *  *  *  *  *  *  * |
         | *  *  *  *  *  *  *  *  * |
         | *  *  *  *  *  *  *  *  * |
         |  *  *  * -a -a -a -a -a -a|
         | *  *  *  *  *  *  *  *  * |
         |-b -b  * -b -b  * -b -b  * |
         | g  g  g  g  g  g  g  g  g |
         | * -b -b  * -b -b  * -b -b |
         | *  *  *  *  *  *  *  *  * |
         |-a -a -a -a -a -a  *  *  * |
         
         Where the repeating cell is of size interiorNy*(interiorNy+1) and is repeated interiorNx times
         */
        
        // Setting of gamma band
        AHat[i*aHatNy + bHatNy*3] = gamma;
        // Setting of lower beta band
        if ((i+1)%bHatNy!=0) {
            AHat[i*aHatNy + 1 + 3*bHatNy] = -1*beta;
        }
        // Setting of upper beta band
        if ((i)%bHatNy!=0) {
            AHat[i*aHatNy + 3*bHatNy-1] = -1*beta;
        }
        // Setting of lower alpha band
        if (i<bHatNy*(bHatNx-1)) {
            AHat[i*aHatNy + 4*bHatNy] = -1*alpha;
        }
        // Setting of upper alpha band
        if (i>bHatNy-1) {
            AHat[i*aHatNy + 2*bHatNy] = -1*alpha;
        }
    }
}

void Poisson2DSolver::GenerateLapackMatrixAHat()
{
    aHatNx = bHatNx*bHatNy;
    aHatNy = (bHatNy*3+1);
    AHat = new double[aHatNx*aHatNy]();
    for (int i=0; i<aHatNx; i++) {
        /* Column-major storage of laplacian banded symmetric matrix for interiorNx=3, interiorNy=3
         | *  *  *  *  *  *  *  *  * |
         | *  *  *  *  *  *  *  *  * |
         | *  *  *  *  *  *  *  *  * |
         |  *  *  * -a -a -a -a -a -a|
         | *  *  *  *  *  *  *  *  * |
         |-b -b  * -b -b  * -b -b  * |
         | g  g  g  g  g  g  g  g  g |
         | * -b -b  * -b -b  * -b -b |
         | *  *  *  *  *  *  *  *  * |
         |-a -a -a -a -a -a  *  *  * |
         
         Where the repeating cell is of size interiorNy*(interiorNy+1) and is repeated interiorNx times
         */
        
        // Setting of gamma band
        AHat[i*aHatNy + bHatNy*2] = gamma;
        // Setting of lower beta band
        if ((i+1)%bHatNy!=0) {
            AHat[i*aHatNy + 1 + 2*bHatNy] = -1*beta;
        }
        // Setting of upper beta band
        if ((i)%bHatNy!=0) {
            AHat[i*aHatNy + 2*bHatNy-1] = -1*beta;
        }
        // Setting of lower alpha band
        if (i<bHatNy*(bHatNx-1)) {
            AHat[i*aHatNy + 3*bHatNy] = -1*alpha;
        }
        // Setting of upper alpha band
        if (i>bHatNy-1) {
            AHat[i*aHatNy + bHatNy] = -1*alpha;
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

void Poisson2DSolver::Updatex(double* xVec) {
    
}

void Poisson2DSolver::SetScalapackMatrixAHat(double* ahat, int ahatnx, int ahatny) {
    delete[] AHat;
    AHat = ahat;
    aHatNx = ahatnx;
    aHatNy = ahatny;
}

void Poisson2DSolver::SolveParallel() {
    
}

void Poisson2DSolver::SolveSerial() {
    int m = bHatNy*bHatNx;      // Number of rows of matrix A
    int n = bHatNy*bHatNx;      // Number of columns of matrix A
    int kl = bHatNy;            // Lower diagonal bandwidth
    int ku = kl;                // Upper diagonal bandwidth
    int lda = 1 + 2*kl + ku;    // Number of rows in compressed matrix
    int* piv = new int[n];      // Pivot data
    int info;
    
    int nrhs = 1;
    
    cblas_dcopy(bHatNy*bHatNx, bHat, 1, xHat, 1);
    F77NAME(dgbtrf) (m, n, kl, ku, AHat, lda, piv, info);
    F77NAME(dgbtrs) ('N', n, kl, ku, nrhs, AHat, lda, piv, xHat, n,info);
}
