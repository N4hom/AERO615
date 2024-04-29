#include <iostream>
#include "Matrix.hpp"
#include "Problem.hpp"
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


	int defN = 40;
	int defM = 120;

	if (argc > 1)
	{
		defN = std::atoi(argv[1]);
	}

	if (argc > 2)
	{
		defM = std::atoi(argv[2]);
	}

	Problem grid = Problem(defN,defM);
	grid.initialize();
	grid.solve();
	grid.save();

	

    auto endTime = std::chrono::high_resolution_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "Run time: " << elapsedTime.count() / 1000.0 << " seconds" << std::endl;

	return 0;
}