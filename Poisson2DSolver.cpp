#include "Poisson2DSolver.h"
#include "cblas.h"
#include <iostream>
#include <iomanip>

Poisson2DSolver::Poisson2DSolver()
{
}

Poisson2DSolver::~Poisson2DSolver()
{
    delete[] A;
    delete[] x;
    delete[] b;
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

void Poisson2DSolver::Initialise(double* xVec, double* bVec, unsigned int bnx, unsigned int bny)
{
    x = xVec;
    b = bVec;
    
    bHatNx = bnx-2;
    bHatNy = bny-2;
}

void Poisson2DSolver::InitialiseMPI()
{
}

void Poisson2DSolver::GenerateScalapackMatrixAHat(double alpha, double beta, double gamma)
{
    aHatNx = bHatNx*bHatNy;
    aHatNy = (bHatNy*2+1);
    A = new double[aHatNx*aHatNy]();
    for (unsigned int i=0; i<bHatNy*bHatNx; i++) {
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
        
        A[i*(bHatNy*2+1) + bHatNy] = gamma;
        if ((i+1)%bHatNy!=0) {
            A[i*(bHatNy*2+1) + 1 + bHatNy] = -1*beta;
        }
        if (i<bHatNy*(bHatNx-1)) {
            A[i*(bHatNy*2+1) + 2*bHatNy] = -1*alpha;
        }
    }
    PrintMatrix2(A,bHatNy,bHatNx,false);
}
