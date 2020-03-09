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
    LidDrivenCavity::Lx = xlen;
    LidDrivenCavity::Ly = ylen;
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
    omega = new double[narr];
    psi = new double[narr];
    
    udy = 2*U/dy;
    narr = Nx*Ny;
        
    cblas_dscal(narr, 0.0, omega, 1);
    cblas_dscal(narr, 0.0, psi, 1);
    
//    cout << Nx;
}

void LidDrivenCavity::Integrate()
{
    SetVorticityBoundaryConditions();
    SetInteriorVorticity();
    UpdateInteriorVorticity();
    SolvePoissonProblem();
}

void LidDrivenCavity::SetVorticityBoundaryConditions()
{
        // Setting top bottom
    for (int i=0; i<Nx; i++) {
        // Setting the top BC
        omega[i+(Ny-1)*Nx] = (psi[i+(Ny-1)*Nx]-psi[i+(Ny-2)*Nx])*2/dy/dy-udy;
        // Setting the bottom BC
        omega[i+(0)*Nx] = (psi[i+(0)*Nx]-psi[i+(1)*Nx])*2/dy/dy;
    }
    
    // Setting left right
    for (int i=0; i<Ny; i++) {
        // Setting the left BC
        omega[i*Nx+0] = (psi[i*Nx+0]-psi[i*Nx+1])*2/dx/dx;
        // Setting the right BC
        omega[i*Nx+(Nx-1)] = (psi[i*Nx+(Nx-1)]-psi[i*Nx+(Nx-2)])*2/dx/dx;
    }
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

void LidDrivenCavity::SetInteriorVorticity()
{
    
}

void LidDrivenCavity::UpdateInteriorVorticity()
{
    
}

void LidDrivenCavity::SolvePoissonProblem()
{
    
}


