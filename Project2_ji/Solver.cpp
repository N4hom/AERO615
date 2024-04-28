// void FlowSolver::correctBoundaryConditions()
// {
//     correctInlet();
//     correctOutlet();
//     correctWall();
// }

// double FlowSolver::correctTimeStep()
// {
//     double minDt = std::numeric_limits<double>::max(); // Start with the largest possible double
    
//     // Loop over each cell to calculate the local time step based on CFL condition
//     for (int i = 0; i < Nci_; ++i) 
//     {
//         for (int j = 0; j < Mci_; ++j) 
//         {

//             // Calculate the denominator of the CFL formula
//             double denom = (lambda_[i][j].s * length_[i][j].s + 
//                             lambda_[i][j].e * length_[i][j].e + 
//                             lambda_[i][j].n * length_[i][j].n + 
//                             lambda_[i][j].w * length_[i][j].w);
//             if (denom != 0) {  // Avoid division by zero
//                 double localDt = CFL_ * (2 * area_[i][j] / denom);
//                 if (localDt < minDt) {
//                     minDt = localDt;
//                 }
//             }
//         }
//     }

//     return minDt;

// }


// void FlowSolver::updateStateProperties() 
// {
//     // Loop over all computational cells
//     for (int i = 0; i < Nc_; ++i) {
//         for (int j = 0; j < Mc_; ++j) {
//             // Update inverse density
//             invRho_[i][j] = 1.0 / q_[0][i][j];  // Assuming q_[0] holds density

//             // Update pressure using the ideal gas law component of the energy equation
//             p_[i][j] = gamma_1_ * (q_[3][i][j] - 0.5 * (pow(q_[1][i][j], 2) + pow(q_[2][i][j], 2)) * invRho_[i][j]);

//             // Update speed of sound
//             c_[i][j] = sqrt(gamma_ * p_[i][j] * invRho_[i][j]);
//         }
//     }
// }


// void FlowSolver::calculateEigen() {
//     for (int i = 0; i < Nci_; ++i) {  // Ensure boundaries are handled, iterate within actual cell boundaries
//         for (int j = 0; j < Mci_; ++j) {
            
            
//             // Used to access the fields with ghost cells at the corresponding i and j
//             int ic = i + 2;
//             int jc = j + 2;

//             std::vector<std::vector<double>>& rhoU = q_[1] ; 
//             std::vector<std::vector<double>>& rhoV = q_[2] ; 

//             // Calculations for each face, using the midpoint formula for averaging between adjacent cells


//             // lambda_[i][j].s = 0.5 * (
//             //     std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_s + rhoV[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].ny_s) +
//             //     std::abs(rhoU[ic    ][jc - 1] * invRho_[ic    ][jc - 1] * n_[i][j].nx_s + rhoV[ic    ][jc - 1] * invRho_[ic    ][jc - 1] * n_[i][j].ny_s))
//             //     + c_[ic    ][jc    ];

//             // // Right face
//             // lambda_[i][j].e = 0.5 * (
//             //     std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_e + rhoV[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].ny_e)  +
//             //     std::abs(rhoU[ic + 1][jc    ] * invRho_[ic + 1][jc    ] * n_[i][j].nx_e + rhoV[ic + 1][jc    ] * invRho_[ic + 1][jc    ] * n_[i][j].ny_e)) 
//             //     + c_[ic    ][jc    ];

//             // // Top face
//             // lambda_[i][j].n = 0.5 * (
//             //     std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_n + rhoV[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].ny_n) + 
//             //     std::abs(rhoU[ic    ][jc + 1] * invRho_[ic    ][jc + 1] * n_[i][j].nx_n + rhoV[ic    ][jc + 1] * invRho_[ic    ][jc + 1] * n_[i][j].ny_n) ) 
//             //     + c_[ic    ][jc    ];

//             // // Left face
//             // lambda_[i][j].w = 0.5 * (
//             //     std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_w + rhoV[ic    ][jc    ] * invRho_[ic     ][jc    ] * n_[i][j].ny_w) + 
//             //     std::abs(rhoU[ic - 1][jc    ] * invRho_[ic - 1][jc    ] * n_[i][j].nx_w + rhoV[ic - 1][jc    ] * invRho_[ic  - 1][jc    ] * n_[i][j].ny_w) ) 
//             //     + c_[ic    ][jc    ];


//             // Bottom face
//             lambda_[i][j].s = 0.5 * (
//                 std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_s + rhoV[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].ny_s) + c_[ic    ][jc    ] +
//                 std::abs(rhoU[ic    ][jc - 1] * invRho_[ic    ][jc - 1] * n_[i][j].nx_s + rhoV[ic    ][jc - 1] * invRho_[ic    ][jc - 1] * n_[i][j].ny_s) + c_[ic    ][jc - 1]
//             );

//             // Right face
//             lambda_[i][j].e = 0.5 * (
//                 std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_e + rhoV[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].ny_e) + c_[ic    ][jc    ] +
//                 std::abs(rhoU[ic + 1][jc    ] * invRho_[ic + 1][jc    ] * n_[i][j].nx_e + rhoV[ic + 1][jc    ] * invRho_[ic + 1][jc    ] * n_[i][j].ny_e) + c_[ic + 1][jc    ]
//             );

//             // Top face
//             lambda_[i][j].n = 0.5 * (
//                 std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_n + rhoV[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].ny_n) + c_[ic    ][jc    ] +
//                 std::abs(rhoU[ic    ][jc + 1] * invRho_[ic    ][jc + 1] * n_[i][j].nx_n + rhoV[ic    ][jc + 1] * invRho_[ic    ][jc + 1] * n_[i][j].ny_n) + c_[ic    ][jc + 1]
//             );

//             // Left face
//             lambda_[i][j].w = 0.5 * (
//                 std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_w + rhoV[ic    ][jc    ] * invRho_[ic     ][jc    ] * n_[i][j].ny_w) + c_[ic    ][jc    ] +
//                 std::abs(rhoU[ic - 1][jc    ] * invRho_[ic - 1][jc    ] * n_[i][j].nx_w + rhoV[ic - 1][jc    ] * invRho_[ic  - 1][jc    ] * n_[i][j].ny_w) + c_[ic - 1][jc    ]
//             );
//         }
//     }
    
//     if (debug_)
//     {
//         printFaceValues(lambda_);
//         /* code */
//     }
// }

// void FlowSolver::runRungeKutta() 
// {
    
//     const std::vector<double> alpha = {0.25, 0.333, 0.5 , 1.0}; 
//     dt_  = correctTimeStep(); // Compute the minimum time step based on CFL condition

//     std::cout << "dt " << dt_ << std::endl;

//     //dt_ = 1e-3;
//     q0_ = q_; // Make a copy of the original state vector

//     for (int stage = 0; stage < 4; ++stage) 
//     {
//         // The most cache-friendly order for the loops would be iterating i 
//         // (over Nci_ + 4), then j (over Mci_ + 4), and finally k (over 4). 
//         // This order ensures that data are accessed that is contiguous in memory, 
//         // which minimizes cache misses and can significantly improve performance, 
//         // especially in large-scale simulations.

//         // Update the state vector for each cell
//         for (int i = 1; i < Nci_ ; ++i) { 
//             for (int j = 0; j < Mci_ ; ++j) {
//                 for (int k = 0; k < 4; ++k) {  // Loop over components

//                     int ic = i + 2;
//                     int jc = j + 2;

//                     q_[k][ic][jc] = q0_[k][ic][jc] - alpha[stage] * dt_ / area_[i][j] * (R_[k][i][j] - D_[k][i][j]);
//                 }
//             }
//         }

//         updateStateProperties(); // Such as 1/rho, pressure, speed of sound
       
//         correctBoundaryConditions();  // Update the BCs if they are dependent on the intermediate states. BCs depend on updated data
//         computeFluxes();   // Update the fluxes based on current state vector
//         computeResiduals(); // Calculate residuals based on the updated fluxes
        
//         // update the state properties if needed
//     }


// }


// void FlowSolver::computeDissipation()
// {
//     calculateEigen();

//     // Loop over each computational cell
//     for (int i = 0; i < Nci_; ++i) {
//         for (int j = 0; j < Mci_; ++j) {
//             // Calculate gradient measures

//             int ic = i + 2;
//             int jc = j + 2;

//             // sCsi and sEta are the switches term in each cell. They are evaluated with a second order central difference scheme
//             double sEta = std::abs(p_[ic + 1][jc    ] - 2 * p_[ic][jc] + p_[ic - 1][jc    ])  / (p_[ic + 1][jc    ] + 2 * p_[ic][jc] + p_[ic - 1][jc    ]);
//             double sCsi = std::abs(p_[ic    ][jc + 1] - 2 * p_[ic][jc] + p_[ic    ][jc - 1])  / (p_[ic + 1][jc    ] + 2 * p_[ic][jc] + p_[ic - 1][jc    ]);

//             // sCsi and sEta evaluated at the east/west (sEta) and north/south (sCsi)
//             // sCsi_i+1/2,j -----> s2_[i][j].s
//             // sCsi_i-1/2,j -----> s2_[i][j].n
//             // sEta_i,j+1/2 -----> s2_[i][j].e
//             // sEta_i,j-1/2 -----> s2_[i][j].w
           
//             s2_[i][j].s = 0.5 * nu2_ * (sCsi + std::abs(p_[ic    ][jc    ] - 2 * p_[ic    ][jc - 1] + p_[ic    ][jc - 2]) / (p_[ic + 1][jc    ] + 2 * p_[ic][jc] + p_[ic - 1][jc    ]));
//             s2_[i][j].n = 0.5 * nu2_ * (sCsi + std::abs(p_[ic    ][jc + 2] - 2 * p_[ic    ][jc + 1] + p_[ic    ][jc    ]) / (p_[ic + 1][jc    ] + 2 * p_[ic][jc] + p_[ic - 1][jc    ]));
//             s2_[i][j].e = 0.5 * nu2_ * (sEta + std::abs(p_[ic + 2][jc    ] - 2 * p_[ic + 1][jc    ] + p_[ic    ][jc    ]) / (p_[ic + 1][jc    ] + 2 * p_[ic][jc] + p_[ic - 1][jc    ]));
//             s2_[i][j].w = 0.5 * nu2_ * (sEta + std::abs(p_[ic    ][jc    ] - 2 * p_[ic - 1][jc    ] + p_[ic - 2][jc    ]) / (p_[ic + 1][jc    ] + 2 * p_[ic][jc] + p_[ic - 1][jc    ]));


//             // Second order dissipation terms
//             s4_[i][j].e = std::max(0.0, nu4_ - s2_[i][j].e);
//             s4_[i][j].w = std::max(0.0, nu4_ - s2_[i][j].w);
//             s4_[i][j].n = std::max(0.0, nu4_ - s2_[i][j].n);  
//             s4_[i][j].s = std::max(0.0, nu4_ - s2_[i][j].s);  
//             // s4_[i][j].n = 0;
//             // s4_[i][j].s = 0;
//             // Compute dissipation terms for each variable
//             for (int k = 0; k < 4; ++k) 
//             {

//                 // (s2.e - s2.w + s2.s - s2.n) -
//                 // (s4.e - s4.w + s2.s - s2.n)

//                 double D2x = s2_[i][j].e * length_[i][j].e * lambda_[i][j].e * (q_[k][ic + 1][jc    ] - q_[k][ic    ][jc    ]) -
//                              s2_[i][j].w * length_[i][j].w * lambda_[i][j].w * (q_[k][ic    ][jc    ] - q_[k][ic - 1][jc    ]);
                
//                 double D2y = (s2_[i][j].n * length_[i][j].n * lambda_[i][j].n * (q_[k][ic    ][jc + 1] - q_[k][ic    ][jc    ])) -
//                              (s2_[i][j].s * length_[i][j].s * lambda_[i][j].s * (q_[k][ic    ][jc    ] - q_[k][ic    ][jc - 1]));
                
//                 double D4x = (s4_[i][j].e * length_[i][j].e * lambda_[i][j].e * (q_[k][ic + 2][jc    ] - 3 * q_[k][ic + 1][jc    ] + 3 * q_[k][ic    ][jc    ] - q_[k][ic - 1][jc    ])) -
//                              (s4_[i][j].w * length_[i][j].w * lambda_[i][j].w * (q_[k][ic + 1][jc    ] - 3 * q_[k][ic    ][jc    ] + 3 * q_[k][ic - 1][jc    ] - q_[k][ic - 2][jc    ]));

//                 double D4y = (s4_[i][j].n * length_[i][j].n * lambda_[i][j].n * (q_[k][ic    ][jc + 2] - 3 * q_[k][ic    ][jc + 1] + 3 * q_[k][ic    ][jc    ] - q_[k][ic    ][jc - 1])) -
//                              (s4_[i][j].s * length_[i][j].s * lambda_[i][j].s * (q_[k][ic    ][jc + 1] - 3 * q_[k][ic    ][jc    ] + 3 * q_[k][ic    ][jc - 1] - q_[k][ic    ][jc - 2]));
                

//                 D_[k][i][j] = (D2x + D2y) - (D4x + D4y);
//             }
//         }
//     }
// }



// void FlowSolver::solve(int iterations, int writeInterval, int verboseInterval) {
//     std::vector<std::vector<std::vector<double>>> dq_max(4, std::vector<std::vector<double>>(Nci_, std::vector<double>(Mci_)));

//     for (int n = 1; n <= iterations; ++n) 
//     {
//         std::cout << "--------------------------------" << std::endl;
//         std::cout << "Iteration  " << n << std::endl;
//         // Calculate dissipation for each cell
//         computeDissipation();

//         if (debug_)
//         {
//                 printVector2D(D_[0] , "rho dissipation ");
//                 printVector2D(D_[1] , "rhoU dissipation ");
//                 printVector2D(D_[2] , "rhoV dissipation ");
//                 printVector2D(D_[3] , "rhoE dissipation ");
//             /* code */
//         }

//         computeResiduals();

//         if (debug_)
//         {
//             printVector2D(R_[0] , "rho dissipation ");
//             printVector2D(R_[1] , "rhoU dissipation ");
//             printVector2D(R_[2] , "rhoV dissipation ");
//             printVector2D(R_[3] , "rhoE dissipation ");
//             /* code */
//         }

//         // Update q using Runge-Kutta 4 steps
//         runRungeKutta();

//         // Calculate maximum residual changes for convergence check
//         for (int k = 0; k < 4; ++k) {
//             for (int i = 0; i < Nci_; ++i) {
//                 for (int j = 0; j < Mci_; ++j) {

//                     int ic = i + 2;
//                     int jc = j + 2;

//                     dq_max[k][i][j] = std::abs(q_[k][ic][jc] - q0_[k][ic][jc]);
//                 }
//             }
//         }


//         // Find global max residuals for each variable
//         std::vector<double> globalMax(4, 0);
//         for (int k = 0; k < 4; ++k) {
//             for (auto& row : dq_max[k]) {
//                 double localMax = *std::max_element(row.begin(), row.end());
//                 if (localMax > globalMax[k]) {
//                     globalMax[k] = localMax;
//                 }
//             }
//         }

//         // Output data and logging
//         if (n % writeInterval == 0) {
//             writeData("STRUCTURED_GRID", n);  // Assuming writeData is correctly implemented
//         }

//         if (n % verboseInterval == 0) {
//             std::cout << "Iteration: " << n;
//             for (auto& max : globalMax) {
//                 std::cout << ", Max Residual: " << max;
//             }
//             std::cout << std::endl;

//             // Assuming file output is setup to write residuals
//             std::ofstream residualFile("residuals.csv", std::ios::app);
//             residualFile << n << ", " << globalMax[0] << ", " << globalMax[1] << ", " << globalMax[2] << ", " << globalMax[3] << "\n";
//         }


//         std::cout << n << ", rho: " << globalMax[0] << ", rhoU" << globalMax[1] << ", rhoV " << globalMax[2] << ", rhoE" << globalMax[3] << "\n";
//                 std::cout << "--------------------------------------------------------------------" << std::endl;

//     }
// }


// void FlowSolver::writeData(const std::string& format, int timestep) {
//     // Construct file name based on the timestep
//     std::string filename = "output_" + std::to_string(timestep) + ".vtk";
//     std::ofstream vtkFile(filename);

//     if (!vtkFile.is_open()) {
//         std::cerr << "Failed to open file for writing VTK data: " << filename << std::endl;
//         return;
//     }

//     // Write VTK header and structured grid format
//     vtkFile << "# vtk DataFile Version 3.0\n";
//     vtkFile << "Structured Grid Example\n";
//     vtkFile << "ASCII\n";
//     vtkFile << "DATASET STRUCTURED_GRID\n";
//     vtkFile << "DIMENSIONS " << Mci_  << " " << Nci_  << " 1\n";
//     vtkFile << "POINTS " << (Mci_) * (Nci_) << " FLOAT\n";

//     // Output grid points assuming a unit spacing between grid points
//     for (int i = 0; i < Nci_ ; i++) {
//         for (int j = 0; j < Mci_ ; j++) {
//             vtkFile << mesh_.cell()[i][j].x  << " " << mesh_.cell()[i][j].y << " 0\n"; // z-coordinate is 0 for 2D grid
//         }
//     }

//     // Write data at points or cells
//     vtkFile << "POINT_DATA " << (Mci_ ) * (Nci_ ) << "\n";
//     vtkFile << "SCALARS pressure float 1\n";
//     vtkFile << "LOOKUP_TABLE default\n";

//     // Output pressure data
//     for (int i = 0; i < Nci_ ; i++) {
//         for (int j = 0; j < Mci_ ; j++) {

//             int ic = i + 2;
//             int jc = j + 2;
//             vtkFile << std::setprecision(6) << p_[ic][jc] << "\n";
//         }
//     }

//     vtkFile << "SCALARS density float 1\n";
//     vtkFile << "LOOKUP_TABLE default\n";
//     for (int i = 0; i < Nci_ ; i++) {
//         for (int j = 0; j < Mci_ ; j++) {
//             int ic = i + 2;
//             int jc = j + 2;
//             vtkFile << std::setprecision(6) << q_[0][ic][jc] << "\n";
//         }
//     }

//     vtkFile << "SCALARS u-velocity float 1\n";
//     vtkFile << "LOOKUP_TABLE default\n";
//     for (int i = 0; i < Nci_ ; i++) {
//         for (int j = 0; j < Mci_ ; j++) {
//             int ic = i + 2;
//             int jc = j + 2;
//             vtkFile << std::setprecision(6) << q_[1][ic][jc]/q_[0][ic][jc] << "\n";
//         }
//     }

//     vtkFile << "SCALARS v-velocity float 1\n";
//     vtkFile << "LOOKUP_TABLE default\n";
//     for (int i = 0; i < Nci_ ; i++) {
//         for (int j = 0; j < Mci_ ; j++) {
//             int ic = i + 2;
//             int jc = j + 2;
//             vtkFile << std::setprecision(6) << q_[2][ic][jc]/q_[0][ic][jc] << "\n";
//         }
//     }

//     vtkFile.close();
//     std::cout << "Data written to " << filename << std::endl;
// }