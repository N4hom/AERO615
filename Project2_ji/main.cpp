#include <iostream>
#include "Mesh.hpp"
#include "Solver.hpp"

int main(int argc, char const *argv[])
{
	Mesh mesh("grids/coarseMeshLeTe.vtk");
	FlowSolver euler(mesh);
	euler.solve(1,1,1);
	return 0;
}