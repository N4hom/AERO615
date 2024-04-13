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




	return 0;
}