#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{
    // Setting up the command line inputs
    po::options_description opts("Solves the lid-driven cavity problem using the given algorithm.");
	opts.add_options()
	    ("help", "Prints help messege.")
		("Lx", po::value<double>()->default_value(1.0), "Set length of the domain in the x-direction. Default = 1.0")
        ("Ly", po::value<double>()->default_value(1.0), "Set length of the domain in the y-direction. Default = 1.0")
        ("Nx", po::value<int>()->default_value(161), "Set number of grid points in the x-direction. Default = 161")
        ("Ny", po::value<int>()->default_value(161), "Set number of grid points in the y-direction. Default = 161")
        ("Px", po::value<int>()->default_value(0), "Set number of partitions in the x-direction. Default = 0")
        ("Py", po::value<int>()->default_value(0), "Set number of partitions in the y-direction. Default = 0")
        ("dt", po::value<double>()->default_value(0.0001), "Set timestep. Default = 1.0")
        ("T", po::value<double>()->default_value(1), "Set final time. Default = 1E-5")
        ("Re", po::value<double>()->default_value(100), "Set Reynolds number. Default = 100");
    
    // Telling boost to parse command-line arguments using list of possible options
    // then generate a map of the options and values
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);
    
    // Checking if user used the "--help" option and print the usage
    if (vm.count("help")) {
        cout << "Solves the lid-driven cavity problem for the given inputs." << endl;
        cout << opts << endl;
        return 0;
    }
    
    // Extracting the options and saving it into variables
    const double Lx = vm["Lx"].as<double>();
    const double Ly = vm["Ly"].as<double>();
    const double Nx = vm["Nx"].as<int>();
    const double Ny = vm["Ny"].as<int>();
    const double Px = vm["Px"].as<int>();
    const double Py = vm["Py"].as<int>();
    const double dt = vm["dt"].as<double>();
    const double T = vm["T"].as<double>();
    const double Re = vm["Re"].as<double>();

    // Displaying the chosen parameters
    cout << "Chosen length of domain in x-direction = " << Lx << endl;
    cout << "Chosen length of domain in y-direction = " << Ly << endl;
    cout << "Chosen number of grid points in x-direction = " << Nx << endl;
    cout << "Chosen number of grid points in y-direction = " << Ny << endl;
    cout << "Chosen number of partitions in x-direction = " << Px << endl;
    cout << "Chosen number of partitions in y-direction = " << Py << endl;
    cout << "Chosen timestep size = " << dt << endl;
    cout << "Chosen final time = " << T << endl;
    cout << "Chosen Reynolds number = " << Re << endl;
    cout << endl;
    
    // Caclulating the grid spacing
    const double dx = Lx / ((double) Nx - 1);
    const double dy = Ly / ((double) Ny - 1);
    
    // Displaying grid spacing
    cout << "The grid space in the x-direction = " << dx << endl;
    cout << "The grid space in the y-direction = " << dy << endl;
    cout << endl;
    
    // Checking the chosen timestep meets the stability restriction
    if (dt >= (Re*dx*dy/4)) {
        cout << "Re*dx*dy/4 = " << Re*dx*dy/4;
        cout << "The time step chosen is too large. Please choose a time step such that dt < Re*dt*dx/4." << endl;
        return 1;
    }
    
    
    // Checking partitions match MPI......
    
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();

    // Configuring Solver
    solver->SetDomainSize(Lx,Ly);
    solver->SetGridSize(Nx,Ny);
    solver->SetTimeStep(dt);
    solver->SetFinalTime(T);
    solver->SetReynoldsNumber(Re);
    solver->SetGridSpacing(dx,dy);
    
    solver->Initialise();
    solver->PrintOmegaMatrix();
    // Run the solver
    // solver->Integrate();
    //
    return 0;
}
