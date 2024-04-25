#include <vector>
#include "Mesh.hpp"

template<class Type>
using Field = std::vector<std::vector<Type>>;


class FlowSolver {
private:
    // Variables
    Mesh& mesh_;
    int icmax_, jcmax_;
    int imax_, jmax_;  // REAL + GHOST CELLS
    int Nc_, Mc_;
    Field<double>& area_;
    Field<FaceNormal>& n_;
    Field<FaceLength>& length_;
    Field<FaceLength>& dx_;
    Field<FaceLength>& dy_;

    double gamma_ = 1.4;
    double gammaM1_ = 0.4;
    double Minf_ = 0.3;
    double pInf_ = 101325;
    double pRatio_ = 0.99;
    double rhoInf_ = 1;
    double cInf_ = sqrt(gamma_ * pInf_ /rhoInf_);
    // std::vector<std::vector<double>> dx; // grid spacing in x-direction
    // std::vector<std::vector<double>> dy; // grid spacing in y-direction
    // std::vector<std::vector<double>> A; // cell area
    // std::vector<std::vector<double>> l; // face length
    // std::vector<std::vector<std::vector<double>>> n; // normal unit vector to face
    
    Field<std::vector<double>> q_; // state vector
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
    Nc_(mesh.Nc()), 
    Mc_(mesh.Mc()), 
    icmax_(mesh.Nc() - 1), 
    jcmax_(mesh.Mc() - 1),
    imax_(icmax_ + 4), 
    jmax_(jcmax_ + 4),
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
        q_.resize(Nc_ + 4, std::vector<std::vector<double>>(Mc_ + 4, std::vector<double>(4)));
        q0_.resize(Nc_ + 4, std::vector<std::vector<double>>(Mc_ + 4, std::vector<double>(4)));
        f_.resize(Nc_ + 4, std::vector<std::vector<double>>(Mc_ + 4, std::vector<double>(4)));
        g_.resize(Nc_ + 4, std::vector<std::vector<double>>(Mc_ + 4, std::vector<double>(4)));
        R_.resize(Nc_, std::vector<std::vector<double>>(Mc_, std::vector<double>(4)));
        D_.resize(Nc_, std::vector<std::vector<double>>(Mc_, std::vector<double>(4)));
        p_.resize(Nc_ + 4, std::vector<double>(Mc_ + 4));
        c_.resize(Nc_ + 4, std::vector<double>(Mc_ + 4));
        invRho_.resize(Nc_ + 4, std::vector<double>(Mc_ + 4));
        lambda_.resize(Nc_, std::vector<double>(Mc_));
        s2_.resize(Nc_, std::vector<double>(Mc_));
        s4_.resize(Nc_, std::vector<double>(Mc_));
    }

    // Method for memory allocation
    void allocateMemory() {
        // Implement memory allocation logic here
    }

    void initializeStateVector(double rhoInf, double uInf, double vInf, double epsInf, double gammaM1) 
    {
        // Initialize state vector q
        for (int i = 0; i < imax_; ++i) {
            for (int j = 0; j < jmax_; ++j) {
                q_[i][j][0] = rhoInf;
                q_[i][j][1] = rhoInf * uInf;
                q_[i][j][2] = rhoInf * vInf;
                q_[i][j][3] = epsInf;
            }
        }

        // Calculate 1 / rho
        for (int i = 0; i < imax_ ; ++i) {
            for (int j = 0; j < jmax_ ; ++j) {
                invRho_[i][j] = 1.0 / q_[i][j][0];
            }
        }

        // Calculate pressure p
        for (int i = 0; i < imax_ ; ++i) {
            for (int j = 0; j < jmax_ ; ++j) {
                p_[i][j] = gammaM1 * (q_[i][j][3] - 0.5 * (q_[i][j][1] * q_[i][j][1] + q_[i][j][2] * q_[i][j][2]) * invRho_[i][j]);
            }
        }

        // Calculate speed of sound c
        for (int i = 0; i < imax_ ; ++i) {
            for (int j = 0; j < jmax_ ; ++j) {
                c_[i][j] = sqrt(gamma * p_[i][j] * invRho_[i][j]);
            }
        }
    }

    // Other methods for initialization, numerical scheme, etc.
};
