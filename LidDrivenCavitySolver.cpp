#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{
    // Setting up the command line inputs
    po::options_description desc("Allowed options");
	desc.add_options()
	    ("help", "produces help messeges")
		("Lx", po::value<double>(), "Set length of the domain in the x-direction")
        ("Ly", po::value<double>(), "Set length of the domain in the y-direction")
        ("Nx", po::value<int>(), "Set number of grid points in the x-direction")
        ("Ny", po::value<int>(), "Set number of grid points in the y-direction")
        ("Px", po::value<int>(), "Set number of partitions in the x-direction")
        ("Py", po::value<int>(), "Set number of partitions in the y-direction")
        ("dt", po::value<double>(), "Set timestep")
        ("T", po::value<double>(), "Set final time")
        ("Re", po::value<double>(), "Set Reynolds number");
        
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        cout << desc << endl;
        return 0;
    }

    if (vm.count("Lx")) {
        cout << "Length of domain in the x-direction was set to " << endl;
    }
    
    if (vm.count("Ly")){
        cout << "Length of domain in the y-direction was set to " << endl;
    }

    if (vm.cout("Nx")){
        cout << "Number of grid points in the x-direction was set to " << endl;
    }
    
    if (vm.count("Ny")){
        cout << "Number of grid points in the y-direction was set to " << endl;
    }

    
    if (vm.count("Px")){
        cout << "Number of partitions in the x-direction was set to " << endl;
    }

 



    // Create a new instance of the LidDrivenCavity class
    //
    // LidDrivenCavity* solver = new LidDrivenCavity();

    // Configure the solver here...
    // ...
    // solver->Initialise();

    // Run the solver
    // solver->Integrate();
    //
    return 0;
}
