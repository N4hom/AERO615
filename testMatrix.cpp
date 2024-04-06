#include <iostream>
#include "Matrix.hpp"

// Read input :
//   - N
//   - Boundary conditions
// Assemble matrix
// Solve
// Output 

int main(int argc, char const *argv[])
{
	Matrix A(5,5,1);
	A.print();

	std::cout << A.rows() << std::endl;

	unsigned int n = 1;
	for (unsigned int i = 0; i < A.rows(); ++i)
	{
		for (int j = 0; j < A.cols(); ++j)
		{
			A(i,j) = n;
			n++;
		}
	}

	A.print() ;


	return 0;
}