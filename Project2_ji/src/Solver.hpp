#ifndef SOLVER_H
#define SOLVER_H 

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <filesystem>
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

inline void printFaceValues(const std::vector<std::vector<Face>>& faceValues)
{
    for (const auto& row : faceValues) {
        for (const auto& elem : row) {
            elem.print();
        }
        std::cout << std::endl;
    }
    std::cout << "\n " << std::endl;
}


struct Vector
{   
    double x;
    double y;

    double mag;

    Vector(double x_coord, double y_coord, double magnitude) : x(x_coord), y(y_coord), mag(magnitude) {};


    void print() const
    {
        std::cout << "( " << x << ", " << y << " )" << " mag: " << mag <<  std::endl;
    }
    
};

class FlowSolver {
private:
    
    // Member data 
    std::string caseName_;

    // Geometry
    double Minf_ ;
    Mesh& mesh_;
    

    int& Nci_;
    int Mci_;
    int Nc_;
    int Mc_;
    int icmax_, jcmax_;
    int imax_, jmax_;  // REAL + GHOST CELLS
    Field<double>& area_;
    Field<FaceNormal>& n_;
    Field<FaceLength>& length_;
    Field<FaceLength>& dx_;
    Field<FaceLength>& dy_;
    Field<Cell>& cell_;
    std::vector<int> bumpIndexBottom_;
    std::vector<int> bumpIndexTop_;


    // Default time step
    double dt_ = 1e-5;

    // double iter_ = 1;
    // double writeInt_ = 1;
    // double logInt_ = 1;


    double CFL_ = 2;
    double nu2_ = 0.1;
    double nu4_ = 0.005;
    double alpha_ = 0;
    double gamma_ = 1.4;
    double gamma_1_ = 0.4;
    double invGammaM1_ = 1/0.4;
    double p0Inf_ = 101325;
    double pRatio_ = 0.99;
    double rhoInf_ = 1.225;
    
    double pInf_ ;
    double cInf_ ;
    double uInf_ ;
    double vInf_ ;
    double epsInf_ ;
    bool isSuperSonic_ = false;
    bool debug_ = false;
    double tolerance_ = 1e-10;

    
    Field<std::vector<double>> q_; // state vector storage
    Field<std::vector<double>> q0_; // state vector at previous time step: used for calculating residuals
    Field<std::vector<double>> f_; // x-component of the fluxes
    Field<std::vector<double>> g_; // y component of the fluxes
    Field<std::vector<double>> R_; // residuals matrix
    Field<std::vector<double>> D_; // dissipation matrix

    // Primitive variables
    Field<double> p_; // pressure
    Field<double> c_; // speed of sound
    
    Field<double> invRho_; // Defined to avoid divisions as it's used multiple time to calculate velocities and speed of sound
    Field<Face> s2_; // second order switch for artificial viscosity 
    Field<Face> lambda_; // eigenvalue
    Field<Face> s4_; // fourth order switch for artificial viscosity
    std::vector<Vector> forceBottom_;
    std::vector<Vector> forceTop_;
    Vector forceBottomIntegral_;
    Vector forceTopIntegral_;

    
    
public:
    
    // Constructor: mesh is passed as a parameter to retrieve mesh info
    FlowSolver(Mesh& mesh, const std::string caseName, const double MachInf);
    
    void initializeStateVector();  // state vector initialization based on provided data
    
    // Boundary conditions
    void correctInlet();  
    void correctOutlet();
    void correctWall();
    void correctBoundaryConditions();
    
    // Flux calculation
    void computeFluxes();

    // Residual and dissipation terms calculation
    void computeResiduals();
    void computeDissipation();

    // Correct time step to satisfy stability
    double correctTimeStep();

    // 4-th order explicit Runge-Kutta scheme
    void runRungeKutta();

    // calculation of the new speed of sound, pressure and inverse of density
    void updateStateProperties();

    // Eigenvalue calculation at each cell face
    void calculateEigen();

    //
    void solve(int iterations, int writeInterval, int residualPrintInterval) ;
    void writeData(int timestep);

    void calcForces();


};


#endif

