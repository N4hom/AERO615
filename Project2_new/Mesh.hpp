#ifndef MESH_H
#define MESH_H


#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "Variable.hpp"

// Define a struct for face length
struct FaceLength {
    double s, e, n, w; // Bottom, Right, Top, Left

    void print() const
    {
        std::cout  << "n = " << n << ", e = " << e << ", s = " << s
                 << ", w = " << w <<  std::endl;

    }
};

struct Face {
    double s, e, n, w; // Bottom, Right, Top, Left

    void print() const
    {
        std::cout  << "n = " << n << ", e = " << e << ", s = " << s
                 << ", w = " << w  << std::endl;

    }
};


// Define a struct for face normal unit vector
struct FaceNormal {
    double nx_s, nx_e, nx_n, nx_w; // Bottom, Right, Top, Left
    double ny_s, ny_e, ny_n, ny_w; // Bottom, Right, Top, Left
};

struct Node {
    double x;
    double y;

    // Constructor
    Node(double x_coord, double y_coord) : x(x_coord), y(y_coord) {}
    friend std::ostream& operator<<(std::ostream& os, const Node& node_) {
        os << "(" << node_.x << ", " << node_.y << ")";
        return os;
    }
};



class Mesh {
private:
    int N_; // Number of nodes in y-direction
    int M_; // Number of nodes in x-direction
    int Nc_; // Number of cells in y-direction
    int Mc_; // Number of cells in x-direction

public:
    // Other members...

    // Method to read grid from file
    void readGrid(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open file." << std::endl;
            return;
        }

        std::string line;
        int inmax, jnmax;
        int nNodes;
        bool foundDataset = false;

        while (std::getline(file, line)) {
            if (line.find("DATASET") != std::string::npos) {
                foundDataset = true;
                break;
            }
        }

        if (!foundDataset) {
            std::cerr << "Error: DATASET line not found." << std::endl;
            return;
        }



        // Read grid dimensions. The dimension show the number of cells in the x direction and then y direction.
        // Since i navigates along rows (y) and j along columns (x), the number of cells in the x direction is the number of columns 
        file >> line >> M_ >> N_;  
        file >> line ;
        std::cout << "Inmax= " << N_ - 1 << ", Jnmax= " << M_ - 1 << std::endl;

        

        // Read number of nodes and verify
        file >> line >> nNodes;
        std::cout << "nNodes : "<< nNodes << std::endl;
        if (N_ * M_ != nNodes) {
            std::cerr << "Error: Number of nodes does not match dimensions." << std::endl;
            return;
        }

        Nc_ = N_ - 1;
        Mc_ = M_ - 1;

        // Resize node_ vector
        node_.resize(N_ , std::vector<Node>(M_ , Node(0.0, 0.0)));

        file >> line;

        for (int i = 0; i < N_; ++i) 
        {
            for (int j = 0; j < M_; ++j) {
                
                file >> node_[i][j].x >> node_[i][j].y; // Swap j and i
                file >> line;
                 std::cout << "Node (" << j << ", " << i << "): x=" << node_[i][j].x << ", y=" << node_[i][j].y << std::endl;
               
            }
        }

        std::cout << "N " << N_ << std::endl;
        std::cout << "M " << M_ << std::endl;
        std::cout << "node_ " << node_[0].size() << std::endl;

        file.close();
    }

    Mesh(const std::string& filename)
    { 
    	readGrid(filename);
    	A_.resize(N_- 1, std::vector<double>(M_ - 1)); // Allocate space once the grid has been read
    	l_.resize(N_- 1, std::vector<FaceLength>(M_ - 1)); // Allocate space once the grid has been read
    	dx_.resize(N_- 1, std::vector<FaceLength>(M_ - 1)); // Allocate space once the grid has been read
    	dy_.resize(N_- 1, std::vector<FaceLength>(M_ - 1)); // Allocate space once the grid has been read
    	n_.resize(N_- 1, std::vector<FaceNormal>(M_ - 1)); // Allocate space once the grid has been read
    	calculateCellProperties();

    	
    }

    // Method to calculate cell properties
    void calculateCellProperties() 
    {
        
        // Variables to store node_ coordinates
        double xA, xB, xC, xD, yA, yB, yC, yD;

        // Loop over cells to calculate properties
        for (int i = 0; i < N_ - 1; ++i) 
        {
            for (int j = 0; j < M_ - 1; ++j) 
            {
                // (A_:left-down node_ CCW to D)
            	//		D_________C
            	//      |		  |
            	//      |		  |
            	//      |		  |
            	//      |_________|
            	//	   A_           B
	            
	            xA = node_[i + 1][j    ].x;
                yA = node_[i + 1][j    ].y;
                xB = node_[i + 1][j + 1].x;
                yB = node_[i + 1][j + 1].y;
                xC = node_[i    ][j + 1].x;
                yC = node_[i    ][j + 1].y;
                xD = node_[i    ][j    ].x;
                yD = node_[i    ][j    ].y;

                // cell area
                A_[i][j] = 0.50 * ((xC - xA) * (yD - yB) - (xD - xB) * (yC - yA));

                // face length
                l_[i][j].s = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA));
                l_[i][j].e = sqrt((xC - xB) * (xC - xB) + (yC - yB) * (yC - yB));
                l_[i][j].n = sqrt((xD - xC) * (xD - xC) + (yD - yC) * (yD - yC));
                l_[i][j].w = sqrt((xA - xD) * (xA - xD) + (yA - yD) * (yA - yD));

                // x component of unit normal vector to the face
                // Assign values to the x and y components of the unit normal vectors
				n_[i][j].nx_s = (yB - yA) / l_[i][j].s; // South face: 
				n_[i][j].ny_s = (xA - xB) / l_[i][j].s; // South face: y component always negative

				n_[i][j].nx_e = (yC - yB) / l_[i][j].e; // East face: x component always positive
				n_[i][j].ny_e = (xB - xC) / l_[i][j].e; // East face

				n_[i][j].nx_n = (yD - yC) / l_[i][j].n; // North face
				n_[i][j].ny_n = (xC - xD) / l_[i][j].n; // North face: y component always positive

				n_[i][j].nx_w = (yA - yD) / l_[i][j].w; // West face: x component always negative
				n_[i][j].ny_w = (xD - xA) / l_[i][j].w; // West face

				// dx of face. It retains the sign for the residual calculation
                dx_[i][j].s = xB - xA;
                dx_[i][j].e = xC - xB;
                dx_[i][j].n = xD - xC;
                dx_[i][j].w = xA - xD;

                // dy of face
                dy_[i][j].s = yB - yA;
                dy_[i][j].e = yC - yB;
                dy_[i][j].n = yD - yC;
                dy_[i][j].w = yA - yD;
	            
            }
        }
    }

    void printNormals()
    {
    	for (int i = 0; i < N_ - 1; ++i) 
        {	
            for (int j = 0; j < M_ - 1; ++j) 
            {
            	std::cout << "i " << i << std::endl;
            	std::cout << "j " << j << std::endl;
            	std::cout << "nx_n " << n_[i][j].nx_n << std::endl;
            	std::cout << "nx_e " << n_[i][j].nx_e << std::endl;
            	std::cout << "nx_s " << n_[i][j].nx_s << std::endl;
            	std::cout << "nx_w " << n_[i][j].nx_w << std::endl;
            	std::cout << "ny_n " << n_[i][j].ny_n << std::endl;
            	std::cout << "ny_e " << n_[i][j].ny_e << std::endl;
            	std::cout << "ny_s " << n_[i][j].ny_s << std::endl;
            	std::cout << "ny_w " << n_[i][j].ny_w << std::endl;
            }
        }
    }

    std::vector<std::vector<Node>>& node(){return node_;}
    std::vector<std::vector<double>>& Area(){return A_;}
    std::vector<std::vector<FaceNormal>>& n(){return n_;}
    std::vector<std::vector<FaceLength>>& length(){return l_;}
    std::vector<std::vector<FaceLength>>& dx(){return dx_;}
    std::vector<std::vector<FaceLength>>& dy(){return dy_;}

    int& Nc(){return Nc_;}
    int& Mc(){return Mc_;}

    

private:
    // Store nodes
    std::vector<std::vector<Node>> node_;

    // Store cell area
    std::vector<std::vector<double>> A_;

    // Store face lengths
    std::vector<std::vector<FaceLength>> l_;

    // Store dx and dy
    std::vector<std::vector<FaceLength>> dx_;
    std::vector<std::vector<FaceLength>> dy_;


    // Store face orientations
    std::vector<std::vector<FaceNormal>> n_;
};


#endif