#include <iostream>
#include <cmath>
#include <iomanip>
#include "cblas.h"

using namespace std;

void printMatrix(double*, int, int, bool);
void printVector(double*, int);

int main(int argc, char **argv)
{
    int Ny = 10;
    int Nx = 15;
    int counter = 1;
    int narr = Nx*Ny;
    
    double* psi = new double[narr];
    double* omega = new double[narr];
    
    for (int i=0; i<Ny; i++) {
        for (int j=0; j<Nx; j++) {
            psi[i*Nx+j] = counter;
            omega[i*Nx+j] = 0;
            counter++;
        }
    }
    
    double dx = 1/((double)Nx-1);
    double dy = 1/((double)Ny-1);
    double udy = 2*1/dy;
    
    cout << dx << " " << dy << endl;
    cout << udy << endl;
    
    // Setting top bottom
    for (int i=0; i<Nx; i++) {
        // Setting the top BC
        omega[i+(Ny-1)*Nx] = (psi[i+(Ny-1)*Nx]-psi[i+(Ny-2)*Nx])*2/dy/dy-udy;
        // Setting the bottom BC
        omega[i+(0)*Nx] = (psi[i+(0)*Nx]-psi[i+(1)*Nx])*2/dy/dy;
    }
    
    // Setting left right
    for (int i=0; i<Ny; i++) {
        // Setting the left BC
        omega[i*Nx+0] = (psi[i*Nx+0]-psi[i*Nx+1])*2/dx/dx;
        // Setting the right BC
        omega[i*Nx+(Nx-1)] = (psi[i*Nx+(Nx-1)]-psi[i*Nx+(Nx-2)])*2/dx/dx;
    }
    
    printMatrix(omega, Ny, Nx, true);
    printMatrix(psi, Ny, Nx, true);
    
//    // Original
//    printMatrix(psi, Ny, Nx, true);
//    printMatrix(omega, Ny, Nx, true);
//
//    // Left
//    cblas_dcopy(narr,psi,Nx,omega,Nx);
//
//    printMatrix(psi, Ny, Nx, true);
//    printMatrix(omega, Ny, Nx, true);
//
//    cblas_daxpy(narr,-1.0,psi+1,Nx,omega,Nx);
//    cblas_dscal()
//
//    printMatrix(psi, Ny, Nx, true);
//    printMatrix(omega, Ny, Nx, true);
    
//    // Bottom
//    cblas_dcopy(narr,psi,Nx,omega,Nx);
//    cblas_daxpy(narr,1.0,psi+1,Nx,omega,Nx);
//
//    printMatrix(psi, Ny, Nx, true);
//    printMatrix(omega, Ny, Nx, true);
    
    cout << "Hello World!" << endl;
    return 0;
}

// Printing Functions

void printMatrix(double* mat, int rows, int cols, bool isRowMajor) {
    for (int i=0; i<rows; i++) {
        for (int j=0; j<cols; j++) {
            if (isRowMajor == false) {
                cout << setw(5) << left << mat[i+j*rows];
            } else {
                cout << setw(5) << left << mat[i*cols+j];
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
