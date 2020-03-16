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
    void SetGridSize(unsigned int nx, unsigned int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void SetGridSpacing(double deltax, double deltay);

    void Initialise();
    void Integrate();
    void GeneratePlots();
    
private:
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
    
    // Discretised grids for omega
    double* omegaInterior = nullptr;
    double* omegaRight = nullptr;
    double* omegaLeft = nullptr;
    double* omegaTop = nullptr;
    double* omegaBottom = nullptr;
    
    // Discretised grids for psi
    double* psiInterior = nullptr;
    double* psiRight = nullptr;
    double* psiLeft = nullptr;
    double* psiTop = nullptr;
    double* psiBottom = nullptr;
    
    // Additional variables
    double udy;
    double alpha;
    double beta;
    double gamma;
    unsigned int narr;
    unsigned int interiorNx;
    unsigned int interiorNy;
    unsigned int interiorNarr;
    double* symmetricBandedLaplacianMatrix = nullptr;
    
    // Private functions for solving vorticity and streamfunction
    void SetVorticityBoundaryConditions();
    void SetInteriorVorticity();
    void UpdateInteriorVorticity();
    
    void SetSymmetricBandedLaplacianMatrix();
};

#endif
