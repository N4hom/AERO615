#include <iostream>
#include "Mesh.hpp"
#include "Solver.hpp"

int main(int argc, char const *argv[])
{
	Mesh mesh("testFine_.vtk");
	FlowSolver euler(mesh);
	euler.solve(10000,5000,5000);
	return 0;
}