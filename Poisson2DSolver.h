/**
    High Performane Computing Coursework
    Poisson2DSolver.h
    Purpose: Defines the Poisson2DSolver class.

    @author Sean Chai
    @version 1.0 23/03/20
*/
#ifndef POISSON2DSOLVER_H
#define POISSON2DSOLVER_H
#pragma once

#include <string>
using namespace std;

/*  For the 2D solver, we assume a problem of [A]{x} = {b}. Since the boundary
    conditions are known, we will reduce this to solving for the interior
    points only. A second order central difference is used for the spatial
    derivatives. The reduced problem becomes [AHat]{xHat} = {bHat}, where
    we apply the known boundary conditions on {b_hat}
 */
class Poisson2DSolver
{
public:
    // Class constructor and destructor
    Poisson2DSolver();
    ~Poisson2DSolver();
    
    // Generation functions
    void GenerateScalapackMatrixAHat();
    void GenerateLapackMatrixAHat();
    
    // Get functions
    double* GetScalapackMatrixAHat();
    int GetScalapackMatrixAHatNx();
    int GetScalapackMatrixAHatNy();
    
    // Set functions
    void SetScalapackMatrixAHat(double* ahat, int ahatnx, int ahatny);
    void SetVariables(int nx, int ny, double alphaVar, double betaVar, double gammaVar);
    void SetVectors(double* xhat, double* bVec);
    void Updatex(double* xVec);
    
    // Prefactor and Solve functions
    void InitialiseScalapack(int px, int py);
    void PrefactorMatrixAHatParallel();
    void SolveParallel();
    void PrefactorMatrixAHatSerial();
    void SolveSerial();
    
private:
    // [A]{x} = {b} matrices and vectors
    double* AHat = nullptr;
    double* xHat = nullptr;
    double* bHat = nullptr;
    
    // Main and superdiagonals of A
    double alpha;
    double beta;
    double gamma;
    
    // Matrix dimensions
    int Nx;
    int Ny;
    int bHatNx;
    int bHatNy;
    int aHatNx;
    int aHatNy;
    
    // MPI variables
    int mpirank;
    int mpisize;
    
    // Scalapack variables
    // BLACS varibles
    int Px;
    int Py;
    int mype;
    int npe;
    int ctx;
    int nrow;
    int ncol;
    int myrow;
    int mycol;
    
    // Solve variables
    int m;
    int n;
    int nb;
    int bwl;
    int bwu;
    int ja;
    int desca[7];
    int lda;
    
    // Solve working arrays
    double* prefactoredAHat = nullptr;
    int* ipiv = nullptr;
    double* af = nullptr;
    int laf;
};

#endif
