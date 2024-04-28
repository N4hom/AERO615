#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "Mesh.hpp"

template<class Type>
using Field = std::vector<std::vector<Type>>;


template<typename T>
void printVector2D(const std::vector<std::vector<T>>& vec2D , const std::string& name="none") {
    std::cout << "Field : " << name << std::endl;
    for (const auto& row : vec2D) {
        for (const auto& elem : row) {
            std::cout << std::scientific << std::setprecision(6) << elem << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\n " << std::endl;
}

template<typename T>
void printStateVector(const std::vector<std::vector<T>>& q) {
    std::vector<std::string> componentNames = {"rho", "rhoU", "rhoV", "rhoE"};
    for (size_t i = 0; i < q.size(); ++i) {
        printVector2D(q[i] , componentNames[i]);
        std::cout << std::endl;
    }
}

template<typename T>
void printFluxes(const std::vector<std::vector<T>>& fg) {
    std::vector<std::string> componentNames = {"rho flux", "rhoU flux", "rhoV flux", "rhoE flux"};
    for (size_t i = 0; i < fg.size(); ++i) {
        printVector2D(fg[i] , componentNames[i]);
        std::cout << std::endl;
    }

}

void printFaceValues(const std::vector<std::vector<Face>>& faceValues)
{
    for (const auto& row : faceValues) {
        for (const auto& elem : row) {
            elem.print();
        }
        std::cout << std::endl;
    }
    std::cout << "\n " << std::endl;
}

class FlowSolver {
private:
    // Variables
    Mesh& mesh_;
    int icmax_, jcmax_;
    int imax_, jmax_;  // REAL + GHOST CELLS
    int& Nci_, Mci_;
    int Nc_, Mc_;
    Field<double>& area_;
    Field<FaceNormal>& n_;
    Field<FaceLength>& length_;
    Field<FaceLength>& dx_;
    Field<FaceLength>& dy_;
    double dt_ = 1e-3;

    double iter_ = 1;
    double writeInt_ = 1;
    double logInt_ = 1;

    double CFL_ = 1;
    double nu2_ = 0.;
    double nu4_ = 0.01;
    double alpha = 0;
    double gamma_ = 1.4;
    double gamma_1_ = 0.4;
    double oneByGammaM1_ = 1/0.4;
    double Minf_ = 0.3;
    double ptInf_ = 101325;
    double pInf_ = ptInf_ / pow(1 + 0.5 * gamma_1_ * Minf_ * Minf_, gamma_*oneByGammaM1_);
    double pRatio_ = 0.99;
    double rhoInf_ = 1.225;
    double cInf_ = sqrt(gamma_ * pInf_ /rhoInf_);
    double uInf_ = Minf_ * cInf_;
    double vInf_ = 0;
    double epsInf_ = pInf_ * oneByGammaM1_ + 0.5 * rhoInf_ * (uInf_ * uInf_ + vInf_ * vInf_);
    bool isSuperSonic_ = false;
    bool debug_ = false;
    
    Field<std::vector<double>> q_; // state vector
    Field<std::vector<double>> qold_; // state vector
    Field<std::vector<double>> q0_; // old state vector
    Field<std::vector<double>> f_; // fluxes
    Field<std::vector<double>> g_; // fluxes
    Field<std::vector<double>> R_; // residuals terms
    Field<std::vector<double>> D_; // dissipation terms
    Field<double> p_; // pressure
    Field<double> c_; // speed of sound
    Field<double> invRho_;
    Field<Face> s2_; // source of 2nd order viscosity
    Field<Face> lambda_; // eigenvalue
    Field<Face> s4_; // source of 4th order viscosity

public:
    // Constructor
    FlowSolver(Mesh& mesh) : 
    mesh_(mesh), 
    Nci_(mesh_.Nc()),  // mesh.Nc() gives the number of real cells
    Mci_(mesh_.Mc()),
    Nc_(mesh_.Nc() + 4),  // Now Nc_ is the number of cells of the computational domain
    Mc_(mesh_.Mc() + 4), 
    icmax_(mesh_.Nc() - 1), 
    jcmax_(mesh_.Mc() - 1),
    imax_(mesh_.Nc() + 3), 
    jmax_(mesh_.Mc() + 3),
    area_(mesh_.Area()),
    n_(mesh_.n()),
    length_(mesh_.length()),
    dx_(mesh_.dx()),
    dy_(mesh_.dy())
    {
        // Resize vectors
        // A.resize(icmax, std::vector<double>(jcmax));
        // l.resize(icmax, std::vector<double>(jcmax));
        // n.resize(icmax, std::vector<std::vector<double>>(jcmax, std::vector<double>(2)));
        //q_.resize(Nci_ + 4, std::vector<std::vector<double>>(Mci_ + 4, std::vector<double>(4)));

        
        std::cout <<"Nc_  " << Nc_ << std::endl;
        std::cout <<"Nci_  " << Nci_ << std::endl;
        std::cout <<"Mci_  " << Mci_ << std::endl;
        std::cout <<"icmax_  " << icmax_ << std::endl;
        std::cout <<"imax_  " << imax_ << std::endl;

        q_.resize(4, std::vector<std::vector<double>>(Nci_ + 4, std::vector<double>(Mci_ + 4)));
        qold_.resize(4, std::vector<std::vector<double>>(Nci_ + 4, std::vector<double>(Mci_ + 4)));
        q0_.resize(4, std::vector<std::vector<double>>(Nci_ + 4, std::vector<double>(Mci_ + 4)));
        f_.resize(4, std::vector<std::vector<double>>(Nci_ + 4, std::vector<double>(Mci_ + 4)));
        g_.resize(4, std::vector<std::vector<double>>(Nci_ + 4, std::vector<double>(Mci_ + 4)));
        R_.resize(4, std::vector<std::vector<double>>(Nci_ , std::vector<double>(Mci_ )));
        D_.resize(4, std::vector<std::vector<double>>(Nci_ , std::vector<double>(Mci_ )));
        p_.resize(Nci_ + 4, std::vector<double>(Mci_ + 4));
        c_.resize(Nci_ + 4, std::vector<double>(Mci_ + 4));
        invRho_.resize(Nci_ + 4, std::vector<double>(Mci_ + 4));
        lambda_.resize(Nci_, std::vector<Face>(Mci_));
        s2_.resize(Nci_, std::vector<Face>(Mci_));
        s4_.resize(Nci_, std::vector<Face>(Mci_));

        


        initializeStateVector();
        if (debug_)
        {
            printStateVector(q_);
            /* code */
        }

        // std::cout << "p_ " << std::endl; 
        // printVector2D(p_ , "p");
        // printStateVector(q_);

        correctInlet();
        if (debug_)
        {
            printStateVector(q_);
            /* code */
        }
        correctOutlet();
        if (debug_)
        {
            printStateVector(q_);
            /* code */
        }
        correctWall();
        
        if (debug_)
        {
            printStateVector(q_);
            /* code */
        }

        updateStateProperties();
        computeFluxes();
        
        if (debug_)
        {
            printVector2D(p_ , "p");
            std::cout << "f fluxes : " << std::endl;
            printFluxes(f_);
            std::cout << "g fluxes : " << std::endl;
            printFluxes(g_);
            /* code */
        }

        //solve(4000 , 1000, 1000);
        

        


        if (Minf_ > 1)
        {
            isSuperSonic_ = true;
        }


    }

    // Method for memory allocation
    void allocateMemory() {
        // Implement memory allocation logic here
    }

    void initializeStateVector() 
    {   
        

        // Loop through until i <= imax or i < N_. Same for j
        // Initialize state vector q
        for (int i = 0; i <= imax_; ++i) 
        {
            for (int j = 0; j <= jmax_; ++j) 
            {
                q_[0][i][j] = rhoInf_;
                q_[1][i][j] = rhoInf_ * uInf_;
                q_[2][i][j] = vInf_;   
                q_[3][i][j] = epsInf_;
            }
        }

        // Calculate 1 / rho
        for (int i = 0; i < Nc_ ; ++i) {
            for (int j = 0; j < Mc_ ; ++j) {
                invRho_[i][j] = 1.0 / q_[0][i][j];
            }
        }

        // Calculate pressure p
        for (int i = 0; i < Nc_ ; ++i) {
            for (int j = 0; j < Mc_ ; ++j) {
                p_[i][j] = gamma_1_ * (q_[3][i][j] - 0.5 * (q_[1][i][j] * q_[1][i][j] + q_[2][i][j] * q_[2][i][j]) * invRho_[i][j]);
            }
        }

        // Calculate speed of sound c
        for (int i = 0; i < Nc_ ; ++i) {
            for (int j = 0; j < Mc_ ; ++j) {
                c_[i][j] = sqrt(gamma_ * p_[i][j] * invRho_[i][j]);
            }
        }
        

    }

    void correctInlet();
    void correctOutlet();
    void correctWall();
    void correctBoundaryConditions();
    void computeFluxes();
    void computeResiduals();
    double correctTimeStep();
    void runRungeKutta();
    void updateStateProperties();
    void calculateEigen();
    void computeDissipation();
    void solve(int iterations, int writeInterval, int verboseInterval) ;
    void writeData(const std::string& format, int timestep);


    // Other methods for initialization, numerical scheme, etc.
};


void FlowSolver::correctInlet() {


    // Loop through all inlet cells
    for (int j = 0; j < Mc_; ++j) {

        int ib = 1;  // Index associated to the inlet

        double u2 = q_[1][ib + 1][j] / q_[0][ib + 1][j];  // u velocity at the inlet
        double v2 = q_[2][ib + 1][j] / q_[0][ib + 1][j];  // v velocity at the inlet
        double c2 = c_[ib + 1][j];  // local speed of sound
        c2 = sqrt(gamma_ * p_[ib + 1][j] * invRho_[ib + 1][j]);
        double Vinf = Minf_ * cInf_;

        // Calculate Riemann invariants (as an example)
        double Riem1 =        Vinf         + 2.0 * cInf_ * oneByGammaM1_;
        double Riem2 = sqrt(u2*u2 + v2*v2) - 2.0 * c2    * oneByGammaM1_;

        // Corrected velocity and sound speed
        double V1 = 0.5 * (Riem1 + Riem2);
        double c1 = 0.25 * gamma_1_ * (Riem1 - Riem2);

        // Calculate pressure at the inlet
        double P1 = ptInf_ * pRatio_/ pow(1.0 + 0.5 * gamma_1_ * (V1 / c1) * (V1 / c1), gamma_ * oneByGammaM1_);
        p_[ib][j] = P1;
        p_[ib - 1][j] = pInf_;

        // Update state vector at the inlet
        q_[0][ib][j] = gamma_ * P1 / (c1 * c1); // Density
        q_[1][ib][j] = q_[0][ib][j] * V1 ;  // Momentum x.  cos(alpha) = 1
        q_[2][ib][j] = 0;  // Momentum y
        q_[3][ib][j] = P1 * oneByGammaM1_ + 0.5 * (pow(q_[1][ib][j], 2) + pow(q_[2][ib][j], 2)) / q_[0][ib][j];  // Energy


        if (isSuperSonic_)
        {
            q_[0][ib][j] = rhoInf_; // Density
            q_[1][ib][j] = rhoInf_ * uInf_ ;  // Momentum x.  cos(alpha) = 1
            q_[2][ib][j] = rhoInf_ * vInf_;  // Momentum y
            q_[3][ib][j] = epsInf_;  // Energy
        }

        q_[0][0][j] = rhoInf_; // Density
        q_[1][0][j] = rhoInf_ * uInf_ ;  // Momentum x.  cos(alpha) = 1
        q_[2][0][j] = rhoInf_ * vInf_;  // Momentum y
        q_[3][0][j] = epsInf_;  // Energy

        // q_[0][0][j] = q_[0][1][j]; // Density
        // q_[1][0][j] = q_[1][1][j];  // Momentum x.  cos(alpha) = 1
        // q_[2][0][j] = q_[2][1][j];  // Momentum y
        // q_[3][0][j] = q_[3][1][j];  // Energy

    }
}



void FlowSolver::correctOutlet()
{
    
    for (int j = 0; j < Mc_; ++j) 
    {
        // I use Jcmax to translate the index by two places
        int Icmax = icmax_ + 2;

        

        std::vector<std::vector<double>>& rho = q_[0];
        std::vector<std::vector<double>>& rhoU = q_[1];
        std::vector<std::vector<double>>& rhoV = q_[2];
        std::vector<std::vector<double>>& rhoE = q_[3];



        rho[Icmax + 1][j]  = 2 * rho[Icmax][j] -   rho[Icmax - 1][j]; // rho(i,0) = rhoInf
        rhoU[Icmax + 1][j] = 2 * rhoU[Icmax][j] - rhoU[Icmax - 1][j]; // rhoU(i,0) = rhoInf * uInf
        rhoV[Icmax + 1][j] = 2 * rhoV[Icmax][j] - rhoV[Icmax - 1][j]; // rhoV(i,0) = rhoInf * vInf
        
        // total energy 
        //rhoE[Icmax + 1][j] = 0.5*pInf_*pRatio_/gamma_1_  + 0.5 * (rhoV[Icmax + 1][j] * rhoV[Icmax + 1][j] + rhoU[Icmax + 1][j] * rhoU[Icmax + 1][j]) / rho[Icmax + 1][j]; // rhoE(i,0) = p/(gamma-1) + 0.5 * rhoInf * uInf * uInf
        rhoE[Icmax + 1][j] = pInf_ * pRatio_ / gamma_1_ + 0.5 * (rhoU[Icmax + 1][j] * rhoU[Icmax + 1][j] + rhoV[Icmax + 1][j] * rhoV[Icmax + 1][j]) / rho[Icmax + 1][j];

        if (isSuperSonic_)
        {
            rho[Icmax + 1 ][j]  = rho[Icmax ][j]; // rho(i,0) = rhoInf
            rhoU[Icmax + 1][j] = rhoU[Icmax ][j]; // rhoU(i,0) = rhoInf * uInf
            rhoV[Icmax + 1][j] = rhoV[Icmax ][j]; // rhoV(i,0) = rhoInf * vInf
            rhoE[Icmax + 1][j] = rhoE[Icmax ][j];
        }

        rho[Icmax + 2][j]  = rho[Icmax + 1][j];
        rhoU[Icmax + 2][j] = rhoU[Icmax + 1][j];
        rhoV[Icmax + 2][j] = rhoV[Icmax + 1][j];
        rhoE[Icmax + 2][j] = rhoE[Icmax + 1][j];
    }
}

void FlowSolver::correctWall()
{

//        j=0                             
//             ________________________________
//             |______|_______|_______|________|
// i = 0       |______|_______|_______|________|     i = icmax_
//             |______|_______|_______|________|
//             |______|_______|_______|________|

//        j=jcmax_                        

    if (debug_)
    {
        std::cout << "Applying wall boundary conditions " << std::endl;
        /* code */
    }
    for (int i = 0; i < Nci_; ++i) 
    {
        // I use Jcmax to translate the index by two places
        int ic = i + 2;
        int Jcmax = jcmax_ + 2;

        

        // q[0] = rho
        // q[1] = rhoU
        // q[2] = rhoV

        std::vector<std::vector<double>>& rho = q_[0];
        std::vector<std::vector<double>>& rhoU = q_[1];
        std::vector<std::vector<double>>& rhoV = q_[2];
        std::vector<std::vector<double>>& rhoE = q_[3];

        // Mirroring the first cells. 
        // data in j = 1 mirrored from j = 2
        rho[ic][1] = rho[ic][2];   // rho(1,j) = rho(2,j)
       
                        // rhoU                 cos(alpha)^2                    sin(alpha)^2                    rhoV            cos(alpha)   sin(alpha)
        rhoU[ic][1] = rhoU[ic][2] * (n_[i][0].ny_s * n_[i][0].ny_s - n_[i][0].nx_s * n_[i][0].nx_s) - 2 *  rhoV[ic][2] * n_[i][0].ny_s * n_[i][0].nx_s;  // Mirroring the u velocity
                                        

                        // rhoV                 sin(alpha)^2                    cos(alpha)^2                    rhoU            cos(alpha)   sin(alpha)
        rhoV[ic][1] = rhoV[ic][2] * (n_[i][0].nx_s * n_[i][0].nx_s - n_[i][0].ny_s * n_[i][0].ny_s) - 2 *  rhoU[ic][2] * n_[i][0].ny_s * n_[i][0].nx_s;  // Mirroring the v velocity
        
        
        
        rhoE[ic][1] = rhoE[ic][2];

        // Mirroring the second layer at the bottom wall
        rho[ic][0] = rho[ic][3];  
       
                        // rhoU                 cos(alpha)^2                    sin(alpha)^2                    rhoV            cos(alpha)   sin(alpha)
        rhoU[ic][0] = rhoU[ic][3] * (n_[i][1].ny_s * n_[i][1].ny_s - n_[i][1].nx_s * n_[i][1].nx_s) - 2 *  rhoV[ic][1] * n_[i][1].ny_s * n_[i][1].nx_s;  // Mirroring the u velocity
        rhoV[ic][0] = rhoV[ic][3] * (n_[i][1].nx_s * n_[i][1].nx_s - n_[i][1].ny_s * n_[i][1].ny_s) - 2 *  rhoU[ic][1] * n_[i][1].ny_s * n_[i][1].nx_s;  // Mirroring the v velocity
                                        
        rhoE[ic][0] = rhoE[ic][3];
        



        // std::cout << "i " << i << std::endl;
        // std::cout << "n_[i][0].nx_n " << n_[i][0].nx_s << std::endl;
        // std::cout << "n_[i][1].ny_n " << n_[i][1].ny_s << std::endl;
        // std::cout << "rhoV[ic][0] " << rhoV[ic][0] << std::endl;
        // std::cout << "rhoV[ic][1] " << rhoV[ic][1] << std::endl;
        // std::cout << "rhoV[ic][2] " << rhoV[ic][2] << std::endl;
        // std::cout << "rhoU[ic][0] " << rhoU[ic][0] << std::endl;
        // std::cout << "rhoU[ic][1] " << rhoU[ic][1] << std::endl;
        // std::cout << "rhoU[ic][1] " << rhoU[ic][2] << std::endl;
        // std::cout << "\n " ;

        ////////////////////////////////////////
        // Mirroring the last cells

        rho[ic][Jcmax + 1] = rho[ic][Jcmax] ; 
        
            // rhoU_                // rhoU 
        rhoU[ic][Jcmax + 1] = rhoU[ic][Jcmax] * (n_[i][jcmax_].ny_n * n_[i][jcmax_].ny_n - n_[i][jcmax_].nx_n * n_[i][jcmax_].nx_n) - 2 *  rhoV[ic][Jcmax] * n_[i][jcmax_].ny_n * n_[i][jcmax_].nx_n;  // Mirroring the u velocity
        rhoV[ic][Jcmax + 1] = rhoV[ic][Jcmax] * (n_[i][jcmax_].nx_n * n_[i][jcmax_].nx_n - n_[i][jcmax_].ny_n * n_[i][jcmax_].ny_n) - 2 *  rhoU[ic][Jcmax] * n_[i][jcmax_].ny_n * n_[i][jcmax_].nx_n;  // Mirroring the v velocity
        

        rhoE[ic][Jcmax + 1] = rhoE[ic][Jcmax];
        
        // Second layer at the bottom wall

        rho[ic][Jcmax + 2] = rho[ic][Jcmax - 1] ; // rho(Icmax + 1,j) = rho(Icmax,j) 
        
            // rhoU_                // rhoU 
        rhoU[ic][Jcmax + 2] = rhoU[ic][Jcmax - 1] * (n_[i][jcmax_ - 1].ny_n * n_[i][jcmax_ - 1].ny_n - n_[i][jcmax_ - 1].nx_n * n_[i][jcmax_ - 1].nx_n) - 2 *  rhoV[i][Jcmax - 1] * n_[i][jcmax_ - 1].ny_n * n_[i][jcmax_ - 1].nx_n;  // Mirroring the u velocity
        rhoV[ic][Jcmax + 2] = rhoV[ic][Jcmax - 1] * (n_[i][jcmax_ - 1].nx_n * n_[i][jcmax_ - 1].nx_n - n_[i][jcmax_ - 1].ny_n * n_[i][jcmax_ - 1].ny_n) - 2 *  rhoU[i][Jcmax - 1] * n_[i][jcmax_ - 1].ny_n * n_[i][jcmax_ - 1].nx_n;  // Mirroring the v velocity
        

        rhoE[ic][Jcmax + 2] = rhoE[ic][Jcmax - 1];


        // std::cout << "i " << i << std::endl;
        // std::cout << "n_[i][jmax].nx_n " << n_[i][jcmax_].nx_n << std::endl;
        // std::cout << "n_[i][jcmax].ny_n " << n_[i][jcmax_].ny_n << std::endl;
        // std::cout << "rhoV[ic][Jcmax    ] " << rhoV[ic][Jcmax    ] << std::endl;
        // std::cout << "rhoV[ic][Jcmax + 1] " << rhoV[ic][Jcmax + 1] << std::endl;
        // std::cout << "rhoV[ic][Jcmax + 2] " << rhoV[ic][Jcmax + 2] << std::endl;
        // std::cout << "rhoU[ic][Jcmax    ] " << rhoU[ic][Jcmax    ] << std::endl;
        // std::cout << "rhoU[ic][Jcmax + 1] " << rhoU[ic][Jcmax + 1] << std::endl;
        // std::cout << "rhoU[ic][Jcmax + 2] " << rhoU[ic][Jcmax + 2] << std::endl;
        // std::cout << "\n " ;
        
        
    }
}


void FlowSolver::computeFluxes() {
    
    // Loop over all cells in the computational domain including ghost cells if needed
    for (int i = 0; i < Nc_; ++i) {
        for (int j = 0; j < Mc_; ++j) {
            
            // Compute fluxes in x-direction (f vector)
            f_[0][i][j] = q_[1][i][j];  // rho * u
            f_[1][i][j] = q_[1][i][j]  *  q_[1][i][j] * invRho_[i][j] + p_[i][j];  // rho * u^2 + p
            f_[2][i][j] = q_[1][i][j]  *  q_[2][i][j] * invRho_[i][j];  // rho * u * v
            f_[3][i][j] = q_[1][i][j]  * (q_[3][i][j] + p_[i][j]) * invRho_[i][j];  // rho * u * (E + p)

            // Compute fluxes in y-direction (g vector)
            g_[0][i][j] = q_[2][i][j];  // rho * v
            g_[1][i][j] = q_[2][i][j]  *  q_[1][i][j] * invRho_[i][j];  // rho * u * v
            g_[2][i][j] = q_[2][i][j]  *  q_[2][i][j] * invRho_[i][j] + p_[i][j];  // rho * v^2 + p
            g_[3][i][j] = q_[2][i][j]  * (q_[3][i][j] + p_[i][j]) * invRho_[i][j];  // rho * v * (E + p)
        }
    }
}

void FlowSolver::computeResiduals() {
    
    // Loop over the components of the state vector
    for (int k = 0; k < 4; ++k) 
    {
        // Loop over all real cells, not including ghost cells
        for (int i = 0; i < Nci_; ++i) 
        { // Adjust index for 0-based and exclude ghost cells
            for (int j = 0; j < Mci_; ++j) 
            {
                int ic = i + 2;
                int jc = j + 2;

                double fW = 0.5 * (f_[k][ic - 1][jc    ] + f_[k][ic][jc]);
                double fE = 0.5 * (f_[k][ic + 1][jc    ] + f_[k][ic][jc]);
                double fN = 0.5 * (f_[k][ic    ][jc + 1] + f_[k][ic][jc]);
                double fS = 0.5 * (f_[k][ic    ][jc - 1] + f_[k][ic][jc]);

                double gW = 0.5 * (g_[k][ic - 1][jc    ] + g_[k][ic][jc]);
                double gE = 0.5 * (g_[k][ic + 1][jc    ] + g_[k][ic][jc]);
                double gN = 0.5 * (g_[k][ic    ][jc + 1] + g_[k][ic][jc]);
                double gS = 0.5 * (g_[k][ic    ][jc - 1] + g_[k][ic][jc]);

                FaceLength& dx_ij = dx_[i][j];
                FaceLength& dy_ij = dy_[i][j];

                R_[k][i][j] = (fE * dy_ij.e + fW * dy_ij.w + fN * dy_ij.n + fS * dy_ij.s) - 
                              (gE * dx_ij.e + gW * dx_ij.w + gN * dx_ij.n + gS * dx_ij.s);

                //R_[k][i][j] = (fE * dy_ij.e + fW * dy_ij.w ) ;  // The f fluxes should be multiplied by dy (e.g. fE *dy.e)

                

                
            }
        }
    }

    if (debug_)
    {
        printVector2D(R_[0] , "rho residual");
        printVector2D(R_[1] , "rhoU residual");
        printVector2D(R_[2] , "rhoV residual");
        printVector2D(R_[3] , "rhoE residual");
    }
}

void FlowSolver::correctBoundaryConditions()
{
    correctInlet();
    correctOutlet();
    correctWall();
}

double FlowSolver::correctTimeStep()
{
    double minDt = std::numeric_limits<double>::max(); // Start with the largest possible double
    
    // Loop over each cell to calculate the local time step based on CFL condition
    for (int i = 0; i < Nci_; ++i) 
    {
        for (int j = 0; j < Mci_; ++j) 
        {

            // Calculate the denominator of the CFL formula
            double denom = (lambda_[i][j].s * length_[i][j].s + 
                            lambda_[i][j].e * length_[i][j].e + 
                            lambda_[i][j].n * length_[i][j].n + 
                            lambda_[i][j].w * length_[i][j].w);
            if (denom != 0) {  // Avoid division by zero
                double localDt = CFL_ * (2 * area_[i][j] / denom);
                if (localDt < minDt) {
                    minDt = localDt;
                }
            }
        }
    }

    return minDt;

}


void FlowSolver::updateStateProperties() 
{
    // Loop over all computational cells
    for (int i = 0; i < Nc_; ++i) {
        for (int j = 0; j < Mc_; ++j) {
            // Update inverse density
            invRho_[i][j] = 1.0 / q_[0][i][j];  // Assuming q_[0] holds density

            // Update pressure using the ideal gas law component of the energy equation
            p_[i][j] = gamma_1_ * (q_[3][i][j] - 0.5 * (pow(q_[1][i][j], 2) + pow(q_[2][i][j], 2)) * invRho_[i][j]);

            // Update speed of sound
            c_[i][j] = sqrt(gamma_ * p_[i][j] * invRho_[i][j]);
        }
    }
}


void FlowSolver::calculateEigen() {
    for (int i = 0; i < Nci_; ++i) {  // Ensure boundaries are handled, iterate within actual cell boundaries
        for (int j = 0; j < Mci_; ++j) {
            
            
            // Used to access the fields with ghost cells at the corresponding i and j
            int ic = i + 2;
            int jc = j + 2;

            std::vector<std::vector<double>>& rhoU = q_[1] ; 
            std::vector<std::vector<double>>& rhoV = q_[2] ; 

            // Calculations for each face, using the midpoint formula for averaging between adjacent cells


            // lambda_[i][j].s = 0.5 * (
            //     std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_s + rhoV[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].ny_s) +
            //     std::abs(rhoU[ic    ][jc - 1] * invRho_[ic    ][jc - 1] * n_[i][j].nx_s + rhoV[ic    ][jc - 1] * invRho_[ic    ][jc - 1] * n_[i][j].ny_s))
            //     + c_[ic    ][jc    ];

            // // Right face
            // lambda_[i][j].e = 0.5 * (
            //     std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_e + rhoV[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].ny_e)  +
            //     std::abs(rhoU[ic + 1][jc    ] * invRho_[ic + 1][jc    ] * n_[i][j].nx_e + rhoV[ic + 1][jc    ] * invRho_[ic + 1][jc    ] * n_[i][j].ny_e)) 
            //     + c_[ic    ][jc    ];

            // // Top face
            // lambda_[i][j].n = 0.5 * (
            //     std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_n + rhoV[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].ny_n) + 
            //     std::abs(rhoU[ic    ][jc + 1] * invRho_[ic    ][jc + 1] * n_[i][j].nx_n + rhoV[ic    ][jc + 1] * invRho_[ic    ][jc + 1] * n_[i][j].ny_n) ) 
            //     + c_[ic    ][jc    ];

            // // Left face
            // lambda_[i][j].w = 0.5 * (
            //     std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_w + rhoV[ic    ][jc    ] * invRho_[ic     ][jc    ] * n_[i][j].ny_w) + 
            //     std::abs(rhoU[ic - 1][jc    ] * invRho_[ic - 1][jc    ] * n_[i][j].nx_w + rhoV[ic - 1][jc    ] * invRho_[ic  - 1][jc    ] * n_[i][j].ny_w) ) 
            //     + c_[ic    ][jc    ];


            // Bottom face
            lambda_[i][j].s = 0.5 * (
                std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_s + rhoV[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].ny_s) + c_[ic    ][jc    ] +
                std::abs(rhoU[ic    ][jc - 1] * invRho_[ic    ][jc - 1] * n_[i][j].nx_s + rhoV[ic    ][jc - 1] * invRho_[ic    ][jc - 1] * n_[i][j].ny_s) + c_[ic    ][jc - 1]
            );

            // Right face
            lambda_[i][j].e = 0.5 * (
                std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_e + rhoV[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].ny_e) + c_[ic    ][jc    ] +
                std::abs(rhoU[ic + 1][jc    ] * invRho_[ic + 1][jc    ] * n_[i][j].nx_e + rhoV[ic + 1][jc    ] * invRho_[ic + 1][jc    ] * n_[i][j].ny_e) + c_[ic + 1][jc    ]
            );

            // Top face
            lambda_[i][j].n = 0.5 * (
                std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_n + rhoV[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].ny_n) + c_[ic    ][jc    ] +
                std::abs(rhoU[ic    ][jc + 1] * invRho_[ic    ][jc + 1] * n_[i][j].nx_n + rhoV[ic    ][jc + 1] * invRho_[ic    ][jc + 1] * n_[i][j].ny_n) + c_[ic    ][jc + 1]
            );

            // Left face
            lambda_[i][j].w = 0.5 * (
                std::abs(rhoU[ic    ][jc    ] * invRho_[ic    ][jc    ] * n_[i][j].nx_w + rhoV[ic    ][jc    ] * invRho_[ic     ][jc    ] * n_[i][j].ny_w) + c_[ic    ][jc    ] +
                std::abs(rhoU[ic - 1][jc    ] * invRho_[ic - 1][jc    ] * n_[i][j].nx_w + rhoV[ic - 1][jc    ] * invRho_[ic  - 1][jc    ] * n_[i][j].ny_w) + c_[ic - 1][jc    ]
            );
        }
    }
    
    if (debug_)
    {
        printFaceValues(lambda_);
        /* code */
    }
}

void FlowSolver::runRungeKutta() 
{
    
    const std::vector<double> alpha = {0.25, 0.333, 0.5 , 1.0}; 
    dt_  = correctTimeStep(); // Compute the minimum time step based on CFL condition

    std::cout << "dt " << dt_ << std::endl;

    //dt_ = 1e-3;
    q0_ = q_; // Make a copy of the original state vector

    for (int stage = 0; stage < 4; ++stage) 
    {
        // The most cache-friendly order for the loops would be iterating i 
        // (over Nci_ + 4), then j (over Mci_ + 4), and finally k (over 4). 
        // This order ensures that data are accessed that is contiguous in memory, 
        // which minimizes cache misses and can significantly improve performance, 
        // especially in large-scale simulations.

        // Update the state vector for each cell
        for (int i = 0; i < Nci_ ; ++i) { 
            for (int j = 0; j < Mci_ ; ++j) {
                for (int k = 0; k < 4; ++k) {  // Loop over components

                    int ic = i + 2;
                    int jc = j + 2;

                    q_[k][ic][jc] = q0_[k][ic][jc] - alpha[stage] * dt_ / area_[i][j] * (R_[k][i][j] - D_[k][i][j]);
                }
            }
        }

        updateStateProperties(); // Such as 1/rho, pressure, speed of sound
       
        correctBoundaryConditions();  // Update the BCs if they are dependent on the intermediate states. BCs depend on updated data
        computeFluxes();   // Update the fluxes based on current state vector
        computeResiduals(); // Calculate residuals based on the updated fluxes
        
        // update the state properties if needed
    }


}


void FlowSolver::computeDissipation()
{
    calculateEigen();

    // Loop over each computational cell
    for (int i = 0; i < Nci_; ++i) {
        for (int j = 0; j < Mci_; ++j) {
            // Calculate gradient measures

            int ic = i + 2;
            int jc = j + 2;

            // sCsi and sEta are the switches term in each cell. They are evaluated with a second order central difference scheme
            double sEta = std::abs(p_[ic + 1][jc    ] - 2 * p_[ic][jc] + p_[ic - 1][jc    ])  / (p_[ic + 1][jc    ] + 2 * p_[ic][jc] + p_[ic - 1][jc    ]);
            double sCsi = std::abs(p_[ic    ][jc + 1] - 2 * p_[ic][jc] + p_[ic    ][jc - 1])  / (p_[ic + 1][jc    ] + 2 * p_[ic][jc] + p_[ic - 1][jc    ]);

            // sCsi and sEta evaluated at the east/west (sEta) and north/south (sCsi)
            // sCsi_i+1/2,j -----> s2_[i][j].s
            // sCsi_i-1/2,j -----> s2_[i][j].n
            // sEta_i,j+1/2 -----> s2_[i][j].e
            // sEta_i,j-1/2 -----> s2_[i][j].w
           
            s2_[i][j].s = 0.5 * nu2_ * (sCsi + std::abs(p_[ic    ][jc    ] - 2 * p_[ic    ][jc - 1] + p_[ic    ][jc - 2]) / (p_[ic + 1][jc    ] + 2 * p_[ic][jc] + p_[ic - 1][jc    ]));
            s2_[i][j].n = 0.5 * nu2_ * (sCsi + std::abs(p_[ic    ][jc + 2] - 2 * p_[ic    ][jc + 1] + p_[ic    ][jc    ]) / (p_[ic + 1][jc    ] + 2 * p_[ic][jc] + p_[ic - 1][jc    ]));
            s2_[i][j].e = 0.5 * nu2_ * (sEta + std::abs(p_[ic + 2][jc    ] - 2 * p_[ic + 1][jc    ] + p_[ic    ][jc    ]) / (p_[ic + 1][jc    ] + 2 * p_[ic][jc] + p_[ic - 1][jc    ]));
            s2_[i][j].w = 0.5 * nu2_ * (sEta + std::abs(p_[ic    ][jc    ] - 2 * p_[ic - 1][jc    ] + p_[ic - 2][jc    ]) / (p_[ic + 1][jc    ] + 2 * p_[ic][jc] + p_[ic - 1][jc    ]));


            // Second order dissipation terms
            s4_[i][j].e = std::max(0.0, nu4_ - s2_[i][j].e);
            s4_[i][j].w = std::max(0.0, nu4_ - s2_[i][j].w);
            s4_[i][j].n = std::max(0.0, nu4_ - s2_[i][j].n);  
            s4_[i][j].s = std::max(0.0, nu4_ - s2_[i][j].s);  
            // s4_[i][j].n = 0;
            // s4_[i][j].s = 0;
            // Compute dissipation terms for each variable
            for (int k = 0; k < 4; ++k) 
            {

                // (s2.e - s2.w + s2.s - s2.n) -
                // (s4.e - s4.w + s2.s - s2.n)

                double D2x = s2_[i][j].e * length_[i][j].e * lambda_[i][j].e * (q_[k][ic + 1][jc    ] - q_[k][ic    ][jc    ]) -
                             s2_[i][j].w * length_[i][j].w * lambda_[i][j].w * (q_[k][ic    ][jc    ] - q_[k][ic - 1][jc    ]);
                
                double D2y = (s2_[i][j].n * length_[i][j].n * lambda_[i][j].n * (q_[k][ic    ][jc + 1] - q_[k][ic    ][jc    ])) -
                             (s2_[i][j].s * length_[i][j].s * lambda_[i][j].s * (q_[k][ic    ][jc    ] - q_[k][ic    ][jc - 1]));
                
                double D4x = (s4_[i][j].e * length_[i][j].e * lambda_[i][j].e * (q_[k][ic + 2][jc    ] - 3 * q_[k][ic + 1][jc    ] + 3 * q_[k][ic    ][jc    ] - q_[k][ic - 1][jc    ])) -
                             (s4_[i][j].w * length_[i][j].w * lambda_[i][j].w * (q_[k][ic + 1][jc    ] - 3 * q_[k][ic    ][jc    ] + 3 * q_[k][ic - 1][jc    ] - q_[k][ic - 2][jc    ]));

                double D4y = (s4_[i][j].n * length_[i][j].n * lambda_[i][j].n * (q_[k][ic    ][jc + 2] - 3 * q_[k][ic    ][jc + 1] + 3 * q_[k][ic    ][jc    ] - q_[k][ic    ][jc - 1])) -
                             (s4_[i][j].s * length_[i][j].s * lambda_[i][j].s * (q_[k][ic    ][jc + 1] - 3 * q_[k][ic    ][jc    ] + 3 * q_[k][ic    ][jc - 1] - q_[k][ic    ][jc - 2]));
                

                D_[k][i][j] = (D2x + D2y) - (D4x + D4y);
            }
        }
    }
}



void FlowSolver::solve(int iterations, int writeInterval, int verboseInterval) {
    std::vector<std::vector<std::vector<double>>> dq_max(4, std::vector<std::vector<double>>(Nci_, std::vector<double>(Mci_)));

    for (int n = 1; n <= iterations; ++n) 
    {
        std::cout << "--------------------------------" << std::endl;
        std::cout << "Iteration  " << n << std::endl;
        // Calculate dissipation for each cell
        computeDissipation();

        if (debug_)
        {
                printVector2D(D_[0] , "rho dissipation ");
                printVector2D(D_[1] , "rhoU dissipation ");
                printVector2D(D_[2] , "rhoV dissipation ");
                printVector2D(D_[3] , "rhoE dissipation ");
            /* code */
        }

        computeResiduals();

        if (debug_)
        {
            printVector2D(R_[0] , "rho dissipation ");
            printVector2D(R_[1] , "rhoU dissipation ");
            printVector2D(R_[2] , "rhoV dissipation ");
            printVector2D(R_[3] , "rhoE dissipation ");
            /* code */
        }

        // Update q using Runge-Kutta 4 steps
        runRungeKutta();

        // Calculate maximum residual changes for convergence check
        for (int k = 0; k < 4; ++k) {
            for (int i = 0; i < Nci_; ++i) {
                for (int j = 0; j < Mci_; ++j) {

                    int ic = i + 2;
                    int jc = j + 2;

                    dq_max[k][i][j] = std::abs(q_[k][ic][jc] - q0_[k][ic][jc]);
                }
            }
        }


        // Find global max residuals for each variable
        std::vector<double> globalMax(4, 0);
        for (int k = 0; k < 4; ++k) {
            for (auto& row : dq_max[k]) {
                double localMax = *std::max_element(row.begin(), row.end());
                if (localMax > globalMax[k]) {
                    globalMax[k] = localMax;
                }
            }
        }

        // Output data and logging
        if (n % writeInterval == 0) {
            writeData("STRUCTURED_GRID", n);  // Assuming writeData is correctly implemented
        }

        if (n % verboseInterval == 0) {
            std::cout << "Iteration: " << n;
            for (auto& max : globalMax) {
                std::cout << ", Max Residual: " << max;
            }
            std::cout << std::endl;

            // Assuming file output is setup to write residuals
            std::ofstream residualFile("residuals.csv", std::ios::app);
            residualFile << n << ", " << globalMax[0] << ", " << globalMax[1] << ", " << globalMax[2] << ", " << globalMax[3] << "\n";
        }


        std::cout << n << ", rho: " << globalMax[0] << ", rhoU" << globalMax[1] << ", rhoV " << globalMax[2] << ", rhoE" << globalMax[3] << "\n";
                std::cout << "--------------------------------------------------------------------" << std::endl;

    }
}


void FlowSolver::writeData(const std::string& format, int timestep) {
    // Construct file name based on the timestep
    std::string filename = "output_" + std::to_string(timestep) + ".vtk";
    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) {
        std::cerr << "Failed to open file for writing VTK data: " << filename << std::endl;
        return;
    }

    // Write VTK header and structured grid format
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Structured Grid Example\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_GRID\n";
    vtkFile << "DIMENSIONS " << Mci_  << " " << Nci_  << " 1\n";
    vtkFile << "POINTS " << (Mci_) * (Nci_) << " FLOAT\n";

    // Output grid points assuming a unit spacing between grid points
    for (int i = 0; i < Nci_ ; i++) {
        for (int j = 0; j < Mci_ ; j++) {
            vtkFile << mesh_.cell()[i][j].x  << " " << mesh_.cell()[i][j].y << " 0\n"; // z-coordinate is 0 for 2D grid
        }
    }

    // Write data at points or cells
    vtkFile << "POINT_DATA " << (Mci_ ) * (Nci_ ) << "\n";
    vtkFile << "SCALARS pressure float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";

    // Output pressure data
    for (int i = 0; i < Nci_ ; i++) {
        for (int j = 0; j < Mci_ ; j++) {

            int ic = i + 2;
            int jc = j + 2;
            vtkFile << std::setprecision(6) << p_[ic][jc] << "\n";
        }
    }

    vtkFile << "SCALARS density float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < Nci_ ; i++) {
        for (int j = 0; j < Mci_ ; j++) {
            int ic = i + 2;
            int jc = j + 2;
            vtkFile << std::setprecision(6) << q_[0][ic][jc] << "\n";
        }
    }

    vtkFile << "SCALARS u-velocity float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < Nci_ ; i++) {
        for (int j = 0; j < Mci_ ; j++) {
            int ic = i + 2;
            int jc = j + 2;
            vtkFile << std::setprecision(6) << q_[1][ic][jc]/q_[0][ic][jc] << "\n";
        }
    }

    vtkFile << "SCALARS v-velocity float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < Nci_ ; i++) {
        for (int j = 0; j < Mci_ ; j++) {
            int ic = i + 2;
            int jc = j + 2;
            vtkFile << std::setprecision(6) << q_[2][ic][jc]/q_[0][ic][jc] << "\n";
        }
    }

    vtkFile << "SCALARS Mach float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < Nci_ ; i++) {
        for (int j = 0; j < Mci_ ; j++) {
            int ic = i + 2;
            int jc = j + 2;
            double Mach_ij = (q_[1][ic][jc] * q_[1][ic][jc] + q_[2][ic][jc]*q_[2][ic][jc])/(q_[0][ic][jc] * q_[0][ic][jc]) / c_[ic][jc] ;
            vtkFile << std::setprecision(6) << Mach_ij << "\n";
        }
    }

    vtkFile.close();
    std::cout << "Data written to " << filename << std::endl;
}




