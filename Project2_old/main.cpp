#include <iostream>
#include "Matrix.hpp"
#include "Mesh.hpp"
#include "Problem.hpp"
#include "debug.hpp"
#include <chrono>
#include <fstream>

// Read input :
//   - N
//   - Boundary conditions
// Assemble matrix
// Solve
// Output 


int main(int argc, char const *argv[])
{	
    auto startTime = std::chrono::high_resolution_clock::now();

    Problem problem(5,5, "x" , "y");

    problem.solve();

    auto endTime = std::chrono::high_resolution_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "Run time: " << elapsedTime.count() / 1000.0 << " seconds" << std::endl;


	return 0;
}