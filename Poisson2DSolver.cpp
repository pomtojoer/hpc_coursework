#include "Poisson2DSolver.h"
#include "cblas.h"
#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <cmath>

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
    
    // For parallel solve
    void F77NAME(pdgbtrf)(const int& n, const int& bwl, const int& bwu,
                          double* A, const int& ja, int* desca,
                          int* ipiv, double* af, int& laf,
                          double* work, int& lwork, int& info);

    void F77NAME(pdgbtrs)(const char& trans, const int& n, const int& bwl,
                          const int& bwu, const int& nrhs, double* A,
                          const int& ja, int* desca, int* ipiv,
                          double* b, int& ib, int* descb, double* af, int& laf,
                          double* work, int& lwork, int& info);

    // BLACS declaration
    void Cblacs_pinfo(int*, int*);
    void Cblacs_get(int, int, int*);
    void Cblacs_gridinit(int*, const char*, int, int);
    void Cblacs_gridinfo(int, int*, int*, int*, int*);
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

void Poisson2DSolver::InitialiseScalapack(int px, int py)
{
    Px = px;
    Py = py;
    
    int MPIInitialised;
    // Checking if MPI was initially initalised
    MPI_Initialized(&MPIInitialised);
    if (!MPIInitialised){
        cout << "Error: MPI was not initialised." << endl;
        throw std::exception();
    } else {
        // Get the rank and comm size on each process.
        MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
    }
    
    char order = 'C';  // block cyclic column major processor mapping
    // Initialises the BLACS world communicator (calls MPI_Init if needed)
    Cblacs_pinfo(&mype, &npe);
    // Get the default system context (i.e. MPI_COMM_WORLD)
    Cblacs_get( 0, 0, &ctx );
    // Initialise a process grid of Py rows and Px columns
    Cblacs_gridinit( &ctx, &order, Py, Px );
    // Get info about the grid to verify it is set up
    Cblacs_gridinfo( ctx, &nrow, &ncol, &myrow, &mycol);
}

void Poisson2DSolver::GenerateScalapackMatrixAHat()
{
    AHat = new double[aHatNx*aHatNy];
    
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
    for (int i=1; i<Nx-1; i++) {
        for (int j=1; j<Ny-1; j++) {
            xVec[i*Ny+j] = xHat[(i-1)*bHatNy+(j-1)];
        }
    }
}

void Poisson2DSolver::SetScalapackMatrixAHat(double* ahat, int ahatnx, int ahatny) {
    AHat = ahat;
    aHatNx = ahatnx;
    aHatNy = ahatny;
}

void Poisson2DSolver::SolveParallel() {
    cout << "Proc " << mpirank << "/" << mpisize << " for MPI, proc " << mype << "/" << npe << " for BLACS in position " << "(" << myrow << "," << mycol << ")/(" << nrow << "," << ncol << ") with local matrix " << endl;
    
    int info; // Status value
    const int n = bHatNy*bHatNx; // Total problem size
    const int nb = (int)ceil(n/(Px*Py)); // Blocking size (number of columns per process)
    const int bwl = bHatNy; // Lower bandwidth
    const int bwu = bHatNy; // Upper bandwidth
    const int nrhs = 1; // Number of RHS to solve
    const int ja = 1; // Start offset in matrix (fortran starts at 1)
    const int ib = 0; // Start offset in RHS vector (fortran starts at 1)
    const int la = (1 + 2*bwl + 2*bwu)*nb;
    int laf = (nb+bwu)*(bwl+bwu)+6*(bwl+bwu)*(bwl+2*bwu); // ScaLAPACK documentation
    
    int desca[7];
    desca[0] = 501;     // banded matrix (1-by-P process grid)
    desca[1] = ctx;     // Context
    desca[2] = n;       // Problem size
    desca[3] = nb; // Blocking of matrix
    desca[4] = 0; // Process row/column
    desca[5] = 1+2*bwl+2*bwu; // Local leading dim
    desca[6] = 0; // Reserved
        
    int descb[7]; // Descriptor for RHS
    descb[0] = 502; // Type
    descb[1] = ctx; // Context
    descb[2] = n; // Problem size
    descb[3] = nb; // Blocking of matrix
    descb[4] = 0; // Process row/column
    descb[5] = nb; // Local leading dim
    descb[6] = 0; // Reserved

//    int lwork = -1;
    int lwork = 10;
    double* work;
    int* ipiv = new int [nb]; // Pivoting array
    double* af = new double[laf];
    
    
    F77NAME(pdgbtrf)(n, bwl, bwu, AHat, ja, desca, ipiv, af, laf, work, lwork, info);
//    lwork = work[0];
//    work = new double[lwork];
//
}

void Poisson2DSolver::SolveSerial() {
    const int m = bHatNy*bHatNx;      // Number of rows of matrix A
    const int n = bHatNy*bHatNx;      // Number of columns of matrix A
    const int kl = bHatNy;            // Lower diagonal bandwidth
    const int ku = kl;                // Upper diagonal bandwidth
    const int lda = 1 + 2*kl + ku;    // Number of rows in compressed matrix
    int* piv = new int[n];      // Pivot data
    int info;
    
    int nrhs = 1;
    
    cblas_dcopy(bHatNy*bHatNx, bHat, 1, xHat, 1);
    F77NAME(dgbtrf) (m, n, kl, ku, AHat, lda, piv, info);
    F77NAME(dgbtrs) ('N', n, kl, ku, nrhs, AHat, lda, piv, xHat, n,info);
}
