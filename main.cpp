#include <iostream>
#include "Matrix.hpp"
#include "Problem.cpp"

// Read input :
//   - N
//   - Boundary conditions
// Assemble matrix
// Solve
// Output 

int main(int argc, char const *argv[])
{
	

	
	Problem grid = Problem(5,5);
	grid.initialize();
	
	grid.solve();


	//std::cout << A.size() << endl;
	

	return 0;
}