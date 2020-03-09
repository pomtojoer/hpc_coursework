#ifndef LIDDRIVENCAVITY_H
#define LIDDRIVENCAVITY_H
#pragma once

#include <string>
using namespace std;

class LidDrivenCavity
{
public:
    LidDrivenCavity();
    ~LidDrivenCavity();

    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void SetGridSpacing(double deltax, double deltay);

    void Initialise();
    void Integrate();
    
    void PrintOmegaMatrix();
    // Add any other public functions

private:
    double* v = nullptr;
    double* s = nullptr;
    
    // Given constants
    const double U = 1.0;

    // Defined parameters
    double dt;
    double T;
    unsigned int Nx;
    unsigned int Ny;
    double Lx;
    double Ly;
    double Re;
    
    // Calculated parameters
    double dx;
    double dy;
    
    // Discretised grids
    double* omega = nullptr;
    double* psi = nullptr;
    
    // Additional constants
    double udy;
    int narr;
    
    void SetVorticityBoundaryConditions();
    void SetInteriorVorticity();
    void UpdateInteriorVorticity();
    void SolvePoissonProblem();
};

#endif
