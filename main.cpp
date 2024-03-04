#include <iostream>
#include "Matrix.hpp"
#include "Problem.cpp"
#include <ctime>
// Read input :
//   - N
//   - Boundary conditions
// Assemble matrix
// Solve
// Output 

int main(int argc, char const *argv[])
{	
	clock_t start;
	start = clock();

	Problem grid = Problem(10,10);
	grid.initialize();
	
	grid.solve2();
	grid.save();

	std::cout << "run time: " << (float)start/CLOCKS_PER_SEC << " seconds"  << std::endl;

	return 0;
}