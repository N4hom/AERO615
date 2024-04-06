#include "Mesh.hpp"

Mesh::Mesh(unsigned int N, unsigned int M, std::string filename):
mesh_(N,M)
{
	read("x")
}

Mesh::~Mesh(){}

Mesh::read(std::string filename)
{
	std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file\n";
        return 1;
    }

    std::string line;
    int i = 0;
    while (std::getline(file, line)) {
        // Split the line into values
        std::vector<double> row = split(line, ',');
        // Add the row to the matrix
        mesh_.setRow(i,row);
        ++i;
    }

    // Close the file
    file.close();
}

std::vector<double> split(const std::string &s, char delimiter) {
    std::vector<double> tokens;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(std::stod(token));
    }
    return tokens;
}