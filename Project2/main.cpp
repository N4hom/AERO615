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
	// std::ifstream file("x");
    // if (!file.is_open()) {
    //     std::cerr << "Unable to open file\n";
    //     return 1;
    // }

    // // Vector to store the matrix
    // //std::vector<std::vector<double>> matrix;

    // Matrix<double> A(5,7);

    // // Read the file line by line
    // std::string line;
    // int i = 0;
    // while (std::getline(file, line)) {
    //     // Split the line into values
    //     std::vector<double> row = split(line, ',');
    //     // Add the row to the matrix
    //     A.setRow(i,row);
    //     ++i;
    // }

    // // Close the file
    // file.close();

    // A.print();

    Mesh mesh(4,5, "x" , "y");

    mesh.x().print();
    mesh.y().print();

    mesh.xC().save("xc");
    mesh.yC().save("yc");

    Problem problem(4,5, "x" , "y");

    problem.solve();

    auto endTime = std::chrono::high_resolution_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "Run time: " << elapsedTime.count() / 1000.0 << " seconds" << std::endl;


	return 0;
}