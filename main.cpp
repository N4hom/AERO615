#include <iostream>
#include "Matrix.hpp"
#include "Problem.cpp"
#include <chrono>
// Read input :
//   - N
//   - Boundary conditions
// Assemble matrix
// Solve
// Output 

int main(int argc, char const *argv[])
{	
	auto startTime = std::chrono::high_resolution_clock::now();

	Problem grid = Problem(40,120);
	grid.initialize();
	grid.solve();
	grid.save();

    auto endTime = std::chrono::high_resolution_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "Run time: " << elapsedTime.count() / 1000.0 << " seconds" << std::endl;

	return 0;
}