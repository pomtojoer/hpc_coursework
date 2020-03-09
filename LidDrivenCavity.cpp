#include "LidDrivenCavity.h"
#include "cblas.h"
#include <iostream>


LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
    delete[] omega;
    delete[] psi;
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    Lx = xlen;
    Ly = ylen;
}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    LidDrivenCavity::Nx = nx;
    LidDrivenCavity::Ny = ny;
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
    LidDrivenCavity::dt = deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
    LidDrivenCavity::T = finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
    LidDrivenCavity::Re = re;
}

void LidDrivenCavity::SetGridSpacing(double deltax, double deltay)
{
    LidDrivenCavity::dx = deltax;
    LidDrivenCavity::dy = deltay;
}

void LidDrivenCavity::Initialise()
{
    const int narr = Nx*Ny;
    
    LidDrivenCavity::omega = new double[narr];
    LidDrivenCavity::psi = new double[narr];
    
    cblas_dscal(narr, 0.0, omega, 1);
    cblas_dscal(narr, 0.0, psi, 1);
}

void LidDrivenCavity::Integrate()
{
}


// Remove these
void LidDrivenCavity::PrintOmegaMatrix() {
    for (int x=0; x<Nx; x++) {
        for (int y=0; y<Ny; y++) {
            cout << omega[x*Nx+y];
        }
        cout << endl;
    }
}
