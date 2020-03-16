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
    void Initialise(double* xVec, double* bVec, unsigned int bnx, unsigned int bny);
    void InitialiseMPI();
    
private:
    double* A = nullptr;
    double* x = nullptr;
    double* b = nullptr;
    
    unsigned int bHatNx;
    unsigned int bHatNy;
    unsigned int aHatNx;
    unsigned int aHatNy;
    
    double* topBC = nullptr;
    double* bottomBC = nullptr;
    double* leftBC = nullptr;
    double* rightBC = nullptr;
    
};

#endif
