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
    Poisson2DSolver();
    ~Poisson2DSolver();
    
    void GenerateScalapackMatrixAHat(double alpha, double beta, double gamma);
    void SetBoundaryConditions(double* topbc, double* bottombc, unsigned int bhatnx, double* leftbc, double* rightbc, unsigned int bhatny);
    void ApplyBoundaryConditions();
    
    
private:
    double* AHat = nullptr;
    double* xHat = nullptr;
    double* bHat = nullptr;
    
    unsigned int AHatNx;
    unsigned int AHatNy;
    unsigned int bHatNx;
    unsigned int bHatNy;
    
    double* topBC = nullptr;
    double* bottomBC = nullptr;
    double* leftBC = nullptr;
    double* rightBC = nullptr;
    
};

#endif
