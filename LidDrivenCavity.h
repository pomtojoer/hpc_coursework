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
    void SetPartitionSize(unsigned int px, unsigned int py);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void SetGridSpacing(double deltax, double deltay);
    void SetMPIConfig();

    void Initialise();
    void Integrate();
    void GeneratePlotData();
    
private:
    // Given constants
    const double U = 1.0;

    // Defined parameters
    double dt;
    double T;
    unsigned int Nx;
    unsigned int Ny;
    unsigned int Px;
    unsigned int Py;
    double Lx;
    double Ly;
    double Re;
    
    // Calculated parameters
    double dx;
    double dy;
    
    // Discretised grids for omega
    double* w;
    double* s;
    
    // Additional variables
    unsigned int narr;
    unsigned int interiorNx;
    unsigned int interiorNy;
    unsigned int interiorNarr;
    
    // Variables for MPI
    int MPIRank;
    int MPISize;
    
    // Variables for individual partitions to work on
    int coordArrLen;
    int* iInnerCoords;
    int* jInnerCoords;
    
    // Variables for poisson solver
    double* scalapackMatrix = nullptr;
    int scalapackMatrixNx;
    int scalapackMatrixNy;
    
    // Private functions for solving vorticity and streamfunction
    void SetVorticityBoundaryConditions();
    void SetInteriorVorticity();
    void UpdateInteriorVorticity();
    
};

#endif
