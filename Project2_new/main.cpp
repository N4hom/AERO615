#include <iostream>
#include "Mesh.hpp"
#include "Solver.hpp"

int main(int argc, char const *argv[])
{
	Mesh mesh("testOneBump.vtk");

	FlowSolver euler(mesh);
	return 0;
}