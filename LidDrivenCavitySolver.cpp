#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
#include <mpi.h>
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
        ("Nx", po::value<unsigned int>()->default_value(161), "Set number of grid points in the x-direction. Default = 161")
        ("Ny", po::value<unsigned int>()->default_value(161), "Set number of grid points in the y-direction. Default = 161")
        ("Px", po::value<unsigned int>()->default_value(1), "Set number of partitions in the x-direction. Default = 1")
        ("Py", po::value<unsigned int>()->default_value(1), "Set number of partitions in the y-direction. Default = 1")
        ("dt", po::value<double>()->default_value(0.0001), "Set timestep. Default = 1E-5")
        ("T", po::value<double>()->default_value(1), "Set final time. Default = 1.0")
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
    const double Nx = vm["Nx"].as<unsigned int>();
    const double Ny = vm["Ny"].as<unsigned int>();
    const double Px = vm["Px"].as<unsigned int>();
    const double Py = vm["Py"].as<unsigned int>();
    const double dt = vm["dt"].as<double>();
    const double T = vm["T"].as<double>();
    const double Re = vm["Re"].as<double>();
    
    // Checking partitions match MPI
    int rank = 0;
    int size = 0;

    // Initialise MPI.
    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        cout << "Error: Failed to initialise MPI" << endl;
        return -1;
    }

    // Get the rank and comm size on each process.
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check number of processes is suitable
    if ((Px*Py) != (unsigned int) size) {
        if (rank == 0) {
            cout << "Error: Number of processes must be Nx * Ny" << endl;
        }
        MPI_Finalize();
        return 0;
    }
    
    if (rank==0) {
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
    }
    
    // Performing checks on input variables
    // Checking that the grid length is valid
    if ((Lx <= 0) || (Ly <= 0)) {
        if (rank==0) {
            cout << "Error: The minimum length of the domain must be greater than 0." << endl;
        }
        return -1;
    }
    
    // Checking that the discretisation is valid
    if (Nx < 3 || Ny < 3) {
        if (rank==0) {
            cout << "Error: The minimum number of discretisation in the x or y direction is 3. Please choose a value greater than or equal to this." << endl;
        }
        return -1;
    }
    
    // Checking that time, final time and reynolds are positive
    if ((T <= 0) || (dt <= 0) || (Re <= 0)) {
        if (rank==0) {
            cout << "Error: The minimum values of T, dt and Re must be greater than 0." << endl;
        }
        return -1;
    }
    
    // Caclulating the grid spacing
    const double dx = Lx / ((double) Nx - 1);
    const double dy = Ly / ((double) Ny - 1);
    
    // Displaying grid spacing
    if (rank==0) {
        cout << "The grid space in the x-direction = " << dx << endl;
        cout << "The grid space in the y-direction = " << dy << endl;
        cout << endl;
    }
    
    // Checking the chosen timestep meets the stability restriction
    if (dt >= (Re*dx*dy/4)) {
        if (rank==0) {
            cout << "Re*dx*dy/4 = " << Re*dx*dy/4;
            cout << "Error: The time step chosen is too large. Please choose a time step such that dt < Re*dt*dx/4." << endl;
        }
        return -1;
    }
    
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();
    
    // Configuring Solver
    solver->SetDomainSize(Lx,Ly);
    solver->SetGridSize(Nx,Ny);
    solver->SetPartitionSize(Px,Py);
    solver->SetTimeStep(dt);
    solver->SetFinalTime(T);
    solver->SetReynoldsNumber(Re);
    solver->SetGridSpacing(dx,dy);
    solver->SetMPIConfig();

    // Running Solver
    solver->Initialise();
    solver->Integrate();
    solver->GeneratePlotData();
    
    // Finalise MPI.
    MPI_Finalize();
    
    return 0;
}
