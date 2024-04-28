#include <iostream>
#include "Mesh.hpp"
#include "Solver.hpp"

int main(int argc, char const *argv[])
{
	Mesh mesh("CoarseMesh.vtk");
	FlowSolver euler(mesh);
	euler.solve(100000,100000,100000);
	return 0;
}