#include <vector>
#include "Mesh.hpp"

template<class Type>
using Field = std::vector<std::vector<Type>>;


template<typename T>
void printVector2D(const std::vector<std::vector<T>>& vec2D , const std::string& name="none") {
    std::cout << "Field : " << name << std::endl;
    for (const auto& row : vec2D) {
        for (const auto& elem : row) {
            std::cout << elem << " ";
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

    double alpha = 0;
    double gamma_ = 1.4;
    double gamma_1_ = 0.4;
    double oneByGammaM1_ = 1/0.4;
    double Minf_ = 0.3;
    double pInf_ = 101325;
    double pRatio_ = 0.99;
    double rhoInf_ = 1;
    double cInf_ = sqrt(gamma_ * pInf_ /rhoInf_);
    double uInf_ = Minf_ * cInf_;
    double vInf_ = 0;
    double epsInf_ = pInf_ * oneByGammaM1_ + 0.5 * rhoInf_ * (uInf_ * uInf_ + vInf_ * vInf_);
    
    
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
    Field<double> s2_; // source of 2nd order viscosity
    Field<double> lambda_; // eigenvalue
    Field<double> s4_; // source of 4th order viscosity

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
    jmax_(mesh_.Nc() + 3),
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

        std::cout <<"here " << std::endl;
        std::cout <<"Nc_  " << Nc_ << std::endl;
        std::cout <<"Nci_  " << Nci_ << std::endl;
        std::cout <<"icmax_  " << icmax_ << std::endl;
        std::cout <<"imax_  " << imax_ << std::endl;
        std::cout <<"mesh_.Nc()  " << mesh_.Nc() << std::endl;

        q_.resize(4, std::vector<std::vector<double>>(Nci_ + 4, std::vector<double>(Mci_ + 4)));
        qold_.resize(Nci_ + 4, std::vector<std::vector<double>>(Mci_ + 4, std::vector<double>(4)));
        q0_.resize(Nci_ + 4, std::vector<std::vector<double>>(Mci_ + 4, std::vector<double>(4)));
        f_.resize(Nci_ + 4, std::vector<std::vector<double>>(Mci_ + 4, std::vector<double>(4)));
        g_.resize(Nci_ + 4, std::vector<std::vector<double>>(Mci_ + 4, std::vector<double>(4)));
        R_.resize(Nci_, std::vector<std::vector<double>>(Mci_, std::vector<double>(4)));
        D_.resize(Nci_, std::vector<std::vector<double>>(Mci_, std::vector<double>(4)));
        p_.resize(Nci_ + 4, std::vector<double>(Mci_ + 4));
        c_.resize(Nci_ + 4, std::vector<double>(Mci_ + 4));
        invRho_.resize(Nci_ + 4, std::vector<double>(Mci_ + 4));
        lambda_.resize(Nci_, std::vector<double>(Mci_));
        s2_.resize(Nci_, std::vector<double>(Mci_));
        s4_.resize(Nci_, std::vector<double>(Mci_));

        std::cout <<"here " << std::endl;


        initializeStateVector();

        std::cout << "p_ " << std::endl; 
        printVector2D(p_ , "p");
        printStateVector(q_);

        correctInlet();
        printStateVector(q_);
        correctOutlet();
        printStateVector(q_);
        correctWall();
        printStateVector(q_);

    }

    // Method for memory allocation
    void allocateMemory() {
        // Implement memory allocation logic here
    }

    void initializeStateVector() 
    {   
        std::cout <<"here " << std::endl;
        std::cout <<"Nc_  " << Nc_ << std::endl;
        std::cout <<"Nci_  " << Nci_ << std::endl;
        std::cout <<"icmax_  " << icmax_ << std::endl;
        std::cout <<"imax_  " << imax_ << std::endl;

        // Loop through until i <= imax or i < N_. Same for j
        // Initialize state vector q
        for (int i = 0; i <= imax_; ++i) {
            for (int j = 0; j <= jmax_; ++j) {
                q_[0][i][j] = rhoInf_;
                q_[1][i][j] = rhoInf_ * uInf_;
                q_[2][i][j] = rhoInf_ * uInf_ ;   // to be changed
                q_[3][i][j] = epsInf_;
            }
        }

        // Calculate 1 / rho
        for (int i = 0; i <= imax_ ; ++i) {
            for (int j = 0; j <= jmax_ ; ++j) {
                invRho_[i][j] = 1.0 / q_[0][i][j];
            }
        }

        // Calculate pressure p
        for (int i = 0; i <= imax_ ; ++i) {
            for (int j = 0; j <= jmax_ ; ++j) {
                p_[i][j] = gamma_1_ * (q_[3][i][j] - 0.5 * (q_[1][i][j] * q_[1][i][j] + q_[2][i][j] * q_[2][i][j]) * invRho_[i][j]);
            }
        }

        // Calculate speed of sound c
        for (int i = 0; i <= imax_ ; ++i) {
            for (int j = 0; j <= jmax_ ; ++j) {
                c_[i][j] = sqrt(gamma_ * p_[i][j] * invRho_[i][j]);
            }
        }
        std::cout <<"here " << std::endl;

    }

    void correctInlet();
    void correctOutlet();
    void correctWall();

    // Other methods for initialization, numerical scheme, etc.
};


void FlowSolver::correctInlet()
{   
    std::cout << "Applying inlet boundary conditions " << std::endl;
    // Inlet Boundary:
    // Set the state vector to the corresponding free stream value
    
    for (int i = 0; i < Nc_; ++i) 
    {
        q_[0][i][0] = rhoInf_ ; // rho(i,0) = rhoInf
        q_[1][i][0] = rhoInf_ * uInf_  ; // rhoU(i,0) = rhoInf * uInf
        q_[2][i][0] = rhoInf_ * uInf_ ; // rhoV(i,0) = rhoInf * vInf  (set to uInf only for debugging)
        q_[3][i][0] = epsInf_ ; // rhoE(i,0) = p/(gamma-1) + 0.5 * rhoInf * uInf * uInf


        q_[0][i][1] = q_[0][i][0];
        q_[1][i][1] = q_[1][i][0];
        q_[2][i][1] = q_[2][i][0];
        q_[3][i][1] = q_[3][i][0];
    }
}


void FlowSolver::correctOutlet()
{
    
    std::cout << "Applying outlet boundary conditions " << std::endl;
    for (int i = 0; i < Nc_; ++i) 
    {
        // I use Jcmax to translate the index by two places
        int Jcmax = jcmax_ + 2;

        std::cout << "jcmax_ " << jcmax_ << std::endl;
        std::cout << "Jcmax "  << Jcmax << std::endl;

        std::vector<std::vector<double>>& rho = q_[0];
        std::vector<std::vector<double>>& rhoU = q_[1];
        std::vector<std::vector<double>>& rhoV = q_[2];
        std::vector<std::vector<double>>& rhoE = q_[3];

        rho[i][Jcmax + 1] = rhoInf_ ; // rho(i,0) = rhoInf
        rhoU[i][Jcmax + 1] = rhoInf_ * uInf_   ; // rhoU(i,0) = rhoInf * uInf
        rhoV[i][Jcmax + 1] = rhoInf_ * vInf_ ; // rhoV(i,0) = rhoInf * vInf
        
        // total energy 
        rhoE[i][Jcmax + 1] = pInf_/gamma_1_  + 0.5 * (q_[1][i][Jcmax + 1] * q_[1][i][Jcmax + 1] + q_[2][i][Jcmax + 1] * q_[2][i][Jcmax + 1]) / q_[0][i][Jcmax + 1]; // rhoE(i,0) = p/(gamma-1) + 0.5 * rhoInf * uInf * uInf


        rho[i][Jcmax + 2] = rho[i][Jcmax + 1];
        rhoU[i][Jcmax + 2] = rhoU[i][Jcmax + 1];
        rhoV[i][Jcmax + 2] = rhoV[i][Jcmax + 1];
        rhoE[i][Jcmax + 2] = rhoE[i][Jcmax + 1];
    }
}

void FlowSolver::correctWall()
{

//        i=0                             
//             ________________________________
//             |______|_______|_______|________|
// j = 0       |______|_______|_______|________|     j = jcmax_
//             |______|_______|_______|________|
//             |______|_______|_______|________|

//        i=icmax_                        


    std::cout << "Applying wall boundary conditions " << std::endl;
    for (int j = 0; j < Mc_; ++j) 
    {
        // I use Jcmax to translate the index by two places
        int jc = j + 2;
        int Icmax = icmax_ + 2;

        std::cout << "icmax_ " << icmax_ << std::endl;
        std::cout << "Icmax "  << Icmax << std::endl;

        // q[0] = rho
        // q[1] = rhoU
        // q[2] = rhoV

        std::vector<std::vector<double>>& rho = q_[0];
        std::vector<std::vector<double>>& rhoU = q_[1];
        std::vector<std::vector<double>>& rhoV = q_[2];
        std::vector<std::vector<double>>& rhoE = q_[3];

        // Mirroring the first cells. The data at i = 0 and i = 1 are set using i = 2 and i = 3
        rho[1][j] = rho[2][j];   // rho(1,j) = rho(2,j)
       
                        // rhoU                 cos(alpha)^2                    sin(alpha)^2                    rhoV            cos(alpha)   sin(alpha)
        rhoU[1][j] = rhoU[2][j] * (n_[0][j].ny_n * n_[0][j].ny_n - n_[0][j].nx_n * n_[0][j].nx_n) - 2 *  rhoV[2][j] * n_[0][j].ny_n * n_[0][j].nx_n;  // Mirroring the u velocity
                                        

                        // rhoV                 sin(alpha)^2                    cos(alpha)^2                    rhoU            cos(alpha)   sin(alpha)
        rhoV[1][j] = rhoV[2][j] * (n_[0][j].nx_n * n_[0][j].nx_n - n_[0][j].ny_n * n_[0][j].ny_n) - 2 *  rhoU[2][j] * n_[0][j].ny_n * n_[0][j].nx_n;  // Mirroring the v velocity
        

        rhoE[1][j] = rhoE[2][j];

        // Mirroring the second layer at the top wall
        rho[0][j] = rho[3][j];   // rho(1,j) = rho(2,j)
       
                        // rhoU                 cos(alpha)^2                    sin(alpha)^2                    rhoV            cos(alpha)   sin(alpha)
        rhoU[0][j] = rhoU[3][j] * (n_[1][j].ny_n * n_[1][j].ny_n - n_[1][j].nx_n * n_[1][j].nx_n) - 2 *  rhoV[1][j] * n_[1][j].ny_n * n_[1][j].nx_n;  // Mirroring the u velocity
                                        

                        // rhoV                 cos(alpha)^2                    sin(alpha)^2                    rhoU            cos(alpha)   sin(alpha)
        rhoV[0][j] = rhoV[3][j] * (n_[1][j].nx_n * n_[1][j].nx_n - n_[1][j].ny_n * n_[1][j].ny_n) - 2 *  rhoU[1][j] * n_[1][j].ny_n * n_[1][j].nx_n;  // Mirroring the v velocity
        

        rhoE[0][j] = rhoE[3][j];

        ////////////////////////////////////////
        // Mirroring the last cells

        rho[Icmax + 1][j] = rho[Icmax][j] ; // rho(Icmax + 1,j) = rho(Icmax,j) 
        
            // rhoU_                // rhoU 
        rhoU[Icmax + 1][j] = rhoU[Icmax][j] * (n_[icmax_][j].ny_s * n_[icmax_][j].ny_s - n_[icmax_][j].nx_s * n_[icmax_][j].nx_s) - 2 *  rhoV[Icmax][j] * n_[icmax_][j].ny_s * n_[icmax_][j].nx_s;  // Mirroring the u velocity
        rhoV[Icmax + 1][j] = rhoV[Icmax][j] * (n_[icmax_][j].nx_s * n_[icmax_][j].nx_s - n_[icmax_][j].ny_s * n_[icmax_][j].ny_s) - 2 *  rhoU[Icmax][j] * n_[icmax_][j].ny_s * n_[icmax_][j].nx_s;  // Mirroring the v velocity
        

        rhoE[Icmax + 1][j] = rhoE[Icmax][j];
        
        // Second layer at the bottom wall

        rho[Icmax + 2][j] = rho[Icmax - 1][j] ; // rho(Icmax + 1,j) = rho(Icmax,j) 
        
            // rhoU_                // rhoU 
        rhoU[Icmax + 2][j] = rhoU[Icmax - 1][j] * (n_[icmax_ - 1][j].ny_s * n_[icmax_ - 1][j].ny_s - n_[icmax_ - 1][j].nx_s * n_[icmax_ - 1][j].nx_s) - 2 *  rhoV[Icmax - 1][j] * n_[icmax_ - 1][j].ny_s * n_[icmax_ - 1][j].nx_s;  // Mirroring the u velocity
        rhoV[Icmax + 2][j] = rhoV[Icmax - 1][j] * (n_[icmax_ - 1][j].nx_s * n_[icmax_ - 1][j].nx_s - n_[icmax_ - 1][j].ny_s * n_[icmax_ - 1][j].ny_s) - 2 *  q_[1][Icmax - 1][j] * n_[icmax_ - 1][j].ny_s * n_[icmax_ - 1][j].nx_s;  // Mirroring the v velocity
        

        rhoE[Icmax + 2][j] = rhoE[Icmax - 1][j];
        
        
    }
}