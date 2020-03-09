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

    double dt;
    double T;
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double Re;
    
    double dx;
    double dy;
    
    double* omega = nullptr;
    double* psi = nullptr;
    
    void Update
};

#endif
