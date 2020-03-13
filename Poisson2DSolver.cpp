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
    delete[] bHat;
    
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

void Poisson2DSolver::GenerateScalapackMatrixAHat(double alpha, double beta, double gamma)
{
    AHat = new double[bHatNx*bHatNy*(bHatNy*2+1)]();
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
        
        AHat[i*(bHatNy*2+1) + bHatNy] = gamma;
        if ((i+1)%bHatNy!=0) {
            AHat[i*(bHatNy*2+1) + 1 + bHatNy] = -1*beta;
        }
        if (i<bHatNy*(bHatNx-1)) {
            AHat[i*(bHatNy*2+1) + 2*bHatNy] = -1*alpha;
        }
    }
}

void Poisson2DSolver::SetBoundaryConditions(double* topbc, double* bottombc, unsigned int bhatnx, double* leftbc, double* rightbc, unsigned int bhatny)
{
    topBC = topbc;
    bottomBC = bottombc;
    leftBC = leftbc;
    rightBC = rightbc;
    
    bHatNx = bhatnx;
    bHatNy = bhatny;
}

void Poisson2DSolver::ApplyBoundaryConditions()
{
    
}
