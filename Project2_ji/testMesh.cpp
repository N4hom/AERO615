#include <iostream>
#include "Mesh.hpp"

int main(int argc, char const *argv[])
{
	Mesh mesh("output.vtk");

	mesh.printNormals();
	return 0;
}