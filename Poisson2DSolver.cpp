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
                          double* b, const int& ib, int* descb, double* af, int& laf,
                          double* work, int& lwork, int& info);
    
    // For sending and receiving matrix
    void F77NAME(dgesd2d)(const int& icontxt, const int& m, const int& n,
                          double* A, const int& lda, const int& rdest,
                          int& cdest);

    void F77NAME(dgerv2d)(const int& icontxt, const int& m, const int& n,
                          double* A, const int& lda, const int& rsrc,
                          int& csrc);


    // BLACS declaration
    void Cblacs_pinfo(int*, int*);
    void Cblacs_get(int, int, int*);
    void Cblacs_gridinit(int*, const char*, int, int);
    void Cblacs_gridinfo(int, int*, int*, int*, int*);
    void Cblacs_gridexit(int);
}

Poisson2DSolver::Poisson2DSolver()
{
}

Poisson2DSolver::~Poisson2DSolver()
{
    delete[] AHat;
    delete[] xHat;
    delete[] bHat;
    
    delete[] prefactoredAHat;
    delete[] ipiv;
    delete[] af;
    
    Cblacs_gridexit(ctx);
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
    
    bHat = new double[bHatNx*bHatNy]{};
    xHat = new double[bHatNx*bHatNy]{};
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
    // Initialise a process grid of 1 rows and Px*Py columns
    Cblacs_gridinit( &ctx, &order, 1, Px*Py );
    // Get info about the grid to verify it is set up
    Cblacs_gridinfo( ctx, &nrow, &ncol, &myrow, &mycol);
}

void Poisson2DSolver::GenerateScalapackMatrixAHat()
{
    // Initialise with all zeros
    aHatNx = bHatNx*bHatNy;
    aHatNy = (bHatNy*4+1);
    AHat = new double[aHatNx*aHatNy]{};
    
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
    AHat = new double[aHatNx*aHatNy]{};
    
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

void Poisson2DSolver::PrefactorMatrixAHatParallel() {
    n = bHatNy*bHatNx;  // (global) Total problem size
    nb = max(bHatNy*2+1, (int)ceil((double)n/(double)npe));    // Blocking size (number of columns per process)
    bwl = bHatNy;       // (global) Lower bandwidth
    bwu = bHatNy;       // (global) Upper bandwidth
    ja = 1;             // Start offset in matrix (fortran starts at 1)
    int info;           // Variable to store success
        
    // Filling desca
    desca[0] = 501;             // banded matrix (1-by-P process grid)
    desca[1] = ctx;             // Context
    desca[2] = n;               // (global) Size of matrix being distributed
    desca[3] = nb;              // Blocking of matrix
    desca[4] = 0;               // Process column of first first column of distributed A
    desca[5] = 1+2*bwl+2*bwu;   // Local leading dim
    desca[6] = 0;               // Reserved
    
    // Filling localAhat
    prefactoredAHat = new double[aHatNy*nb]{};
    int offseta = nb*mype*aHatNy;
    for (int i=0; i<nb; i++) {
        for (int j=0; j<aHatNy; j++) {
            if ((i*aHatNy+j+offseta) < (aHatNy*aHatNx)) {
                prefactoredAHat[i*aHatNy+j] = AHat[i*aHatNy+j+offseta];
            }
        }
    }

    // Generating af
    laf = (nb+bwu)*(bwl+bwu)+6*(bwl+bwu)*(bwl+2*bwu);
    af = new double[laf];
    
    // Query for optimal workspace
    int lwork = -1;
    double workopt;
    
    F77NAME(pdgbtrf)(n, bwl, bwu, prefactoredAHat, ja, desca, ipiv, af, laf, &workopt, lwork, info);
    
    // allocating optimal workspace
    lwork = (int)workopt;
    double* work = new double[lwork];
    
    // Generating ipiv. Minimum size of nb+bwu+bwl (bug described in forum http:/icl.cs.utk.edu/lapack-forum/viewtopic.php?f=13&t=2243
    // Original documentation gives nb, but this was found to be buggy with larger array sizes
    ipiv = new int[nb+bwu+bwl];
    // Factorising and solving matrix
    F77NAME(pdgbtrf)(n, bwl, bwu, prefactoredAHat, ja, desca, ipiv, af, laf, work, lwork, info);

    if (info) {
        cout << "Error occurred in PDGBTRF: " << info << endl;
    }
    delete[] work;
}

void Poisson2DSolver::SolveParallel() {
//    cout << "Proc " << mpirank << "/" << mpisize << " for MPI, proc " << mype << "/" << npe << " for BLACS in position " << "(" << myrow << "," << mycol << ")/(" << nrow << "," << ncol << ") with local matrix " << endl;
    
    const int nrhs = 1; // Number of RHS to solve
    double* localbhat;  // Local array b
    const int ib = 1;   // Start offset in RHS vector (fortran starts at 1)
    int descb[7];       // Descriptor for RHS
    int info;           // Status value
    
    // Filling local Bhat
    localbhat = new double[nb]{};
    int offsetb = nb*mype;
    for (int i=0; i<nb; i++) {
        if ((i+offsetb) < (bHatNx*bHatNy)) {
            localbhat[i] = bHat[i+offsetb];
        }
    }
    
    // Filling descb
    descb[0] = 502; // Type
    descb[1] = ctx; // Context
    descb[2] = n;   // Problem size
    descb[3] = nb;  // Blocking of matrix
    descb[4] = 0;   // Process row/column
    descb[5] = nb;  // Local leading dim
    descb[6] = 0;   // Reserved
    
    
    int lwork = (nb+bwu)*(bwl+bwu)+6*(bwl+bwu)*(bwl+2*bwu);
    double* work = new double[lwork];
    F77NAME(pdgbtrs)('N', n, bwl, bwu, nrhs, prefactoredAHat, ja, desca, ipiv, localbhat, ib, descb, af, laf, work, lwork, info);
    
    if (info) {
        cout << "Error occurred in PDGBTRS: " << info << endl;
    }
    
    if (mype > 0) {
        int rdest = 0;
        int cdest = 0;
        F77NAME(dgesd2d)(ctx, nb, 1, localbhat, nb, rdest, cdest);
    } else {
        double* assembledB = new double[nb*npe]{};
        cblas_dcopy(nb, localbhat, 1, assembledB, 1);
        for (int src=1; src<npe; src++) {
            F77NAME(dgerv2d)(ctx, nb, 1, assembledB+src*nb, nb, 0, src);
        }
        cblas_dcopy(bHatNx*bHatNy, assembledB, 1, xHat, 1);
        delete[] assembledB;
    }
    
    delete[] work;
    delete[] localbhat;
}

void Poisson2DSolver::PrefactorMatrixAHatSerial() {
    m = bHatNy*bHatNx;      // Number of rows of matrix A
    n = bHatNy*bHatNx;      // Number of columns of matrix A
    bwl = bHatNy;            // Lower diagonal bandwidth
    bwu = bwl;                // Upper diagonal bandwidth
    lda = 1 + 2*bwl + bwu;    // Number of rows in compressed matrix
    ipiv = new int[n];      // Pivot data
    int info;
    
    F77NAME(dgbtrf) (m, n, bwl, bwu, AHat, lda, ipiv, info);
    if (info) {
        cout << "Error occurred in DGBTRF: " << info << endl;
    }
}

void Poisson2DSolver::SolveSerial() {
    int info;
    int nrhs = 1;
    
    cblas_dcopy(bHatNy*bHatNx, bHat, 1, xHat, 1);
    F77NAME(dgbtrs) ('N', n, bwl, bwu, nrhs, AHat, lda, ipiv, xHat, n, info);
    if (info) {
        cout << "Error occurred in DGBTRS: " << info << endl;
    }
}
