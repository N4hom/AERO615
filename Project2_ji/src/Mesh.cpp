#include "Mesh.hpp"


Mesh::Mesh(const std::string& filename)
{ 
        readGrid(filename);
        A_.resize(N_- 1, std::vector<double>(M_ - 1)); // Allocate space once the grid has been read
        l_.resize(N_- 1, std::vector<FaceLength>(M_ - 1)); 
        dx_.resize(N_- 1, std::vector<FaceLength>(M_ - 1)); 
        dy_.resize(N_- 1, std::vector<FaceLength>(M_ - 1)); 
        n_.resize(N_ - 1, std::vector<FaceNormal>(M_ - 1 ));
        c_.resize(N_- 1, std::vector<Cell>(M_ - 1 , Cell(0.0 , 0.0))); 
        calculateCellProperties();

        findBumps();
            
}

void Mesh::readGrid(const std::string& filename)
{
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open file." << std::endl;
            return;
        }

        std::string line;
        
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
        file >> line >> N_ >> M_;  
        file >> line ;
        std::cout << "Inmax= " << N_ - 1 << ", Jnmax= " << M_ - 1 << std::endl;

        

        // Read number of nodes and verify
        file >> line >> nNodes;
        std::cout << "nNodes : "<< nNodes << std::endl;
        std::cout << "nCells : "<< (N_ - 1) * (M_ - 1) << std::endl;
        if (N_ * M_ != nNodes) {
            std::cerr << "Error: Number of nodes does not match dimensions." << std::endl;
            return;
        }

        Nc_ = N_ - 1;
        Mc_ = M_ - 1;

        // Resize node_ vector
        node_.resize(N_ , std::vector<Node>(M_ , Node(0.0, 0.0)));

        file >> line;


        std::cout << std::fixed << std::setprecision(16);

        {
            for (int j = 0; j < M_; ++j) 
                {
                    // Looping keeping y constant
                    for (int i = 0; i < N_; ++i) {
                        
                        file >> node_[i][j].x >> node_[i][j].y; // Swap j and i
                        file >> line;
                         
                         if (debug_)
                         {
                            std::cout << "Node (" << i << ", " << j << "): x=" << node_[i][j].x << ", y=" << node_[i][j].y << std::endl;
                         }
                       
                    }
                }
        }

        

        file.close();
}


void Mesh::findBumps()
{

    int nCells = c_.size();

    
    for (int i = 0; i < nCells; ++i)
    {
       
        double xBottom = c_[i][0       ].x;
        double xTop = c_[i][M_ - 2].x;

        if (xBottom >= 2 && xBottom <= 3)
        {
            bumpIndexBottom_.push_back(i);           
        }
        
        if (xTop >= 2 && xTop <= 3)
        {
        
            bumpIndexTop_.push_back(i);
        }
    }

    
   
}

void Mesh::calculateCellProperties()
{
        
        // Variables to store node_ coordinates
        double xA, xB, xC, xD, yA, yB, yC, yD;

        // Loop over cells to calculate properties
        for (int i = 0; i < N_ - 1; ++i) 
        {
            for (int j = 0; j < M_ - 1; ++j) 
            {
                // (counter-clockwise)
                //      D_________C
                //      |         |
                //      |         |
                //      |         |
                //      |_________|
                //     A_           B
                
                xA = node_[i    ][j    ].x;
                yA = node_[i    ][j    ].y;
                xB = node_[i + 1][j    ].x;
                yB = node_[i + 1][j    ].y;
                xC = node_[i + 1][j + 1].x;
                yC = node_[i + 1][j + 1].y;
                xD = node_[i    ][j + 1].x;
                yD = node_[i    ][j + 1].y;

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

                // Calculate centroids
                c_[i][j].x = (xA + xB + xC + xD) / 4.0;
                c_[i][j].y = (yA + yB + yC + yD) / 4.0;
                
            }
        }
    }

void Mesh::printNormals()
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

void Mesh::printWallNormals()
{
        
            for (int i = 0; i < N_ - 1; ++i) 
            {
                std::cout << "i " << i << std::endl;
                std::cout << "nx_s " << n_[i][0].nx_s << std::endl;
                std::cout << "ny_s " << n_[i][0].ny_s << std::endl;
            }

            std::cout << "\n ";

            for (int i = 0; i < N_ - 1; ++i) 
            {
                std::cout << "i " << i << std::endl;
                std::cout << "nx_s " << n_[i][1].nx_s << std::endl;
                std::cout << "ny_s " << n_[i][1].ny_s << std::endl;
            }
        
}