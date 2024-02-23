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
	Matrix A(2,2,1);
	A.print();

	
	Problem grid = Problem(5,5);
	grid.initialize();
	grid.x().print();
	grid.y().print();
	grid.alpha().print();
	grid.beta().print();
	grid.gamma().print();


	//std::cout << A.size() << endl;
	

	return 0;
}