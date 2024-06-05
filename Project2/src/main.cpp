#include <iostream>
#include "Mesh.hpp"
#include "Solver.hpp"

int main(int argc, char const *argv[])
{
	// Check if the correct number of arguments are passed
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <VTK file> <case name> <double value>" << std::endl;
        return 1;
    }

    // Parse the command line arguments
    std::string vtkFileName = argv[1];
    std::string caseName = argv[2];
    double MachInf = std::stod(argv[3]); 


	Mesh mesh(vtkFileName);
	FlowSolver euler(mesh , caseName , MachInf);

	//        iterations   output     residuals and forces
	euler.solve(1000000,     100000,      100);
	
	euler.calcForces();
	return 0;
}