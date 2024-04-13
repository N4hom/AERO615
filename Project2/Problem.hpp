#include "Matrix.hpp"
#include "Mesh.hpp"
#include "Variable.hpp"
#include "debug.hpp"
#include <math.h>

class Problem
{

public:
	Mesh mesh_;
	unsigned int Nci_ , Mci_;
	unsigned int Nc_ , Mc_;

	Variable        rho_, rhoU_,rhoV_, rhoE_;  // State vector components
	Matrix<double> U_, V_ , p_, T_ , Umag_;							   // Primitive variables 
	Matrix<double> c_;
	Matrix<double> M_;
	Matrix<double> Riem1_;
	Matrix<double> Riem2_;

	// Matrix<double> Rrho_;
	// Matrix<double> RrhoU_;
	// Matrix<double> RrhoV_;
	// Matrix<double> RrhoE_;

	// Matrix<double> Drho_;
	// Matrix<double> DrhoU_;
	// Matrix<double> DrhoV_;
	// Matrix<double> DrhoE_;

	double gamma_ = 1.4;
	double pInf_  = 101325;;
	double Minf_  = 0.3;
	double alphaInf_ = 0;
	double rhoInf_   = 1;
	double cInf_     = sqrt(gamma_*pInf_/rhoInf_);
	double Uinf_  = Minf_ * cInf_;  

	double deltaCsi(const Matrix<double>& matrix, unsigned int i, unsigned int j);
	double deltaEta(const Matrix<double>& matrix, unsigned int i, unsigned int j);
	double delta2Csi(const Matrix<double>& matrix, unsigned int i, unsigned int j);
	double delta2Eta(const Matrix<double>& matrix, unsigned int i, unsigned int j);
	double delta3Csi(const Matrix<double>& matrix, unsigned int i, unsigned int j);
	double delta3Eta(const Matrix<double>& matrix, unsigned int i, unsigned int j);


	void initialize();
	void initializeFluxes();
	void correctBoundaryConditions();
	void correctInlet();
	void correctWall();
	void correctPeriodic();
	void correctOutlet();
	void correctProperties();
	void R();
	Matrix<double> R(const Matrix<double>& flux_f , const Matrix<double>& flux_g);
	Matrix<double> D(const Matrix<double>& variable);
	void D(const Matrix<double>& flux_f , const Matrix<double>& flux_g);
	void solve();
	void RungeKutta(Variable& variable , int i , int j);
	
	Problem(unsigned int N, unsigned int M, std::string filename_x, std::string filename_y);
	~Problem();
	
};

// Input is the number of nodes
Problem::Problem(unsigned int N, unsigned int M, std::string filename_x, std::string filename_y):
mesh_(N, M, filename_x, filename_y),
Nci_(N - 1),
Mci_(M - 1),
Nc_(Nci_ + 4),
Mc_(Mci_ + 4),
rho_(Nc_ , Mc_ , mesh_),  // N - 1 being the number of internal cells
rhoU_(Nc_ , Mc_, mesh_),
rhoV_(Nc_ , Mc_, mesh_),
rhoE_(Nc_ , Mc_, mesh_),
U_(Nc_ , Mc_),
V_(Nc_ , Mc_),
p_(Nc_ , Mc_),
T_(Nc_ , Mc_),
c_(Nc_ , Mc_),
M_(Nc_ , Mc_),
Umag_(Nc_ , Mc_),
Riem1_(Nc_ , Mc_),
Riem2_(Nc_ , Mc_)
{

	initialize();
	rho_.print();
	rhoU_.print();
	rhoV_.print();
	rhoE_.print();
	
	initializeFluxes();
	

	if(DEBUG)
	{	
		std::cout<< "rho flux f:" << std::endl;
		rho_.flux_f().print();
		std::cout<< "rhoU flux f:" << std::endl;
		rhoU_.flux_f().print();
		std::cout<< "rhoV flux f:" << std::endl;
		rhoV_.flux_f().print();
		std::cout<< "rhoE flux f:" << std::endl;
		rhoE_.flux_f().print();
	}

	
	if (DEBUG)
	{
		std::cout<< "rho flux g:" << std::endl;
		rho_.flux_g().print();
		std::cout<< "rhoU flux g:" << std::endl;
		rhoU_.flux_g().print();
		std::cout<< "rhoV flux g:" << std::endl;
		rhoV_.flux_g().print();
		std::cout<< "rhoE flux g:" << std::endl;
		rhoE_.flux_g().print();
	}

	
}

Problem::~Problem(){}


double Problem::deltaCsi(const Matrix<double>& matrix, unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	double matrixTop = 0.5 * (matrix(ic , jc) + matrix(ic + 1, jc));
	double matrixBottom = 0.5 * (matrix(ic , jc) + matrix(ic - 1, jc));

	return matrixTop - matrixBottom;
}

double Problem::deltaEta(const Matrix<double>& matrix, unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	double matrixRight = 0.5 * (matrix(ic , jc) + matrix(ic , jc + 1));
	double matrixLeft = 0.5 * (matrix(ic , jc) + matrix(ic , jc - 1));

	return matrixRight - matrixLeft;
}

double Problem::delta2Csi(const Matrix<double>& matrix, unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	
	return matrix(ic + 1, jc) - 2*matrix(ic , jc) + matrix(ic - 1, jc);
}

double Problem::delta2Eta(const Matrix<double>& matrix, unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	
	return matrix(ic , jc + 1) - 2*matrix(ic , jc) + matrix(ic , jc + 1);
}

void Problem::initialize()
{

	const unsigned int Nc = mesh_.Nc_;
	const unsigned int Mc = mesh_.Mc_;


	// Initialize internal field
	for (unsigned int i = 0; i < Nc; ++i)
	{
		for (unsigned int j = 0; j < Mc; ++j)
		{
			unsigned ic = i + 2;
			unsigned jc = j + 2;
			rho_(ic , jc) = rhoInf_;
			p_(ic , jc) = pInf_;
			c_(ic , jc) = cInf_;
			U_(ic , jc) = Uinf_;
			rhoU_(ic , jc) = rho_(ic , jc) * U_(ic , jc);
			rhoV_(ic , jc) = rho_(ic , jc) * V_(ic , jc);
			Umag_(ic , jc) = sqrt(U_(ic , jc)*U_(ic , jc) + V_(ic , jc) * V_(ic , jc));
			rhoE_(ic , jc) = p_(ic , jc)/(gamma_ - 1) + 0.5 * Umag_(ic , jc)*Umag_(ic , jc);
		}
	}


	// Initialize inlet boundary
	for (unsigned int i = 0; i < Nc; ++i)
	{
		for (int j = 0 ; j < 3; ++j)
		{
				
			unsigned ic = i + 2; 
			// I don't use jc here because there's no need to since I'm setting the boundary conditions to
			// the initial values. Only the first two columns of cells are touched.

			rho_(ic , j)  = rhoInf_;
			p_(ic , j)    = pInf_;
			c_(ic , j)    = cInf_;
			U_(ic , j )   = Minf_ * cInf_;
			rhoU_(ic , j) = rho_(ic , j) * U_(ic , j);
			rhoV_(ic , j) = rho_(ic , j) * V_(ic , j);
			Umag_(ic , j) = sqrt(U_(ic , j)*U_(ic , j) + V_(ic , j) * V_(ic , j));
			rhoE_(ic , j) = p_(ic , j) / (gamma_ - 1) + 0.5 * Umag_(ic , j)*Umag_(ic , j);
		}
		
	}

	correctWall();
	correctOutlet();
}

void Problem::initializeFluxes()
{
	// f flux component
	Matrix<double>& rhoFlux_f(rho_.flux_f());
	Matrix<double>& rhoUFlux_f(rhoU_.flux_f());
	Matrix<double>& rhoVFlux_f(rhoV_.flux_f());
	Matrix<double>& rhoEFlux_f(rhoE_.flux_f());

	// g flux component
	Matrix<double>& rhoFlux_g(rho_.flux_g());
	Matrix<double>& rhoUFlux_g(rhoU_.flux_g());
	Matrix<double>& rhoVFlux_g(rhoV_.flux_g());
	Matrix<double>& rhoEFlux_g(rhoE_.flux_g());


	// Looping over internal domain to initialize the fluxes in the interior
	for (unsigned int i = 0; i < Nci_; ++i)
	{
		for (unsigned int j = 0; j < Mci_; ++j)
		{
			unsigned ic = i + 2;
			unsigned jc = j + 2;


			rhoFlux_f(ic , jc)  = rhoU_ (ic , jc);
			rhoUFlux_f(ic , jc) = rhoU_(ic , jc) * U_(ic , jc) + p_(ic , jc);
			rhoVFlux_f(ic , jc) = rhoV_(ic , jc) * U_(ic , jc);
			rhoEFlux_f(ic , jc) = rhoE_(ic , jc) * U_(ic , jc) + p_(ic , jc) * U_(ic , jc);
			
			rhoFlux_g(ic , jc)  = rhoV_ (ic , jc);
			rhoUFlux_g(ic , jc) = rhoU_(ic , jc) * V_(ic , jc);
			rhoVFlux_g(ic , jc) = rhoV_(ic , jc) * V_(ic , jc) + p_(ic , jc);
			rhoEFlux_g(ic , jc) = rhoE_(ic , jc) * V_(ic , jc) + p_(ic , jc) * V_(ic , jc);
			
		}
	}

}


void Problem::solve()
{
	rho_.computeResidual();
	rhoU_.computeResidual();
	rhoV_.computeResidual();
	rhoE_.computeResidual();


	for (int i = 0; i < Nci_; ++i)
	{
		for (int j = 0; j < Mci_; ++j)
		{
			RungeKutta(rho_  , i , j);
			RungeKutta(rhoU_ , i , j);
			RungeKutta(rhoV_ , i , j);
			RungeKutta(rhoE_ , i , j);
		}
	}

	correctProperties();
	
	correctWall();
	correctInlet();
	correctOutlet();

	if (DEBUG)
	{
		std::cout << "rho " << std::endl;
		rho_.print();
		std::cout << "rhoU " << std::endl;
		rhoU_.print();
		std::cout << "rhoV " << std::endl;
		rhoV_.print();
		std::cout << "rhoE " << std::endl;
		rhoE_.print();

		std::cout << "c " << std::endl;
		c_.print();
	}
}


void Problem::RungeKutta(Variable& variable, int i , int j )
{
	double alpha1 = 0.25, alpha2 = 0.5 , alpha3 = 1/3, alpha4 = 1;
	double dt = 1e-5;

	Matrix<double>& area = mesh_.area_;

	int ic = i + 2;
	int jc = j + 2;

	double& Aij = area(ic,jc);
	double& Rij =  variable.R()(i , j);
	double  Dij0 = variable.D()(i , j);
	double variable0 = variable(ic , jc);

	// At the next RK step the fluxes must be updated to calculate the new residual. It will be changed only the flux and residual at the cell (ic jc) or (i j)
	
	// Step 1
	variable(ic , jc) = variable(ic , jc) - alpha1 * dt/Aij * (Rij - Dij0);


	// Step 2
	variable.correctFlux_f(U_ , i , j);
	variable.correctFlux_g(V_ , i , j);
	variable(ic , jc) = variable0 - alpha2 * dt/Aij * (Rij - Dij0);


	// Step 3
	variable.correctFlux_f(U_ , i , j);
	variable.correctFlux_g(V_ , i , j);
	variable(ic , jc) = variable0 - alpha3 * dt/Aij * (Rij - Dij0);


	// Step 4
	variable.correctFlux_f(U_ , i , j);
	variable.correctFlux_g(V_ , i , j);
	variable(ic , jc) = variable0 - alpha4 * dt/Aij * (Rij - Dij0);


	


}

// void Problem::initialize()
// {

// 	const unsigned int Nc = mesh_.Nc_;
// 	const unsigned int Mc = mesh_.Mc_;


// 	// Initialize internal field
// 	for (unsigned int i = 0; i < Nc; ++i)
// 	{
// 		for (unsigned int j = 0; j < Mc; ++j)
// 		{
// 			unsigned ic = i + 2;
// 			unsigned jc = j + 2;
// 			rho_(ic , jc) = rhoInf_;
// 			p_(ic , jc) = pInf_;
// 			c_(ic , jc) = cInf_;
// 			U_(ic , jc) = Uinf_;
// 			Umag_(ic , jc) = sqrt(U_(ic , jc)*U_(ic , jc) + V_(ic , jc) * V_(ic , jc));
// 			rhoE_(ic , jc) = p_(ic , jc)/(gamma_ - 1) + 0.5 * Umag_(ic , jc)*Umag_(ic , jc);
// 		}
// 	}


// 	// Initialize inlet boundary
// 	for (unsigned int i = 0; i < Nc; ++i)
// 	{
// 		for (int j = 0 ; j < 3; ++j)
// 		{
				
// 			unsigned ic = i + 2; 
// 			// I don't use jc here because there's no need to since I'm setting the boundary conditions to
// 			// the initial values. Only the first two columns of cells are touched.

// 			rho_(ic , j)  = rhoInf_;
// 			p_(ic , j)    = pInf_;
// 			c_(ic , j)    = cInf_;
// 			U_(ic , j )   = Minf_ * cInf_;
// 			rhoU_(ic , j) = rho_(ic , j) * U_(ic , j);
// 			rhoV_(ic , j) = rho_(ic , j) * V_(ic , j);
// 			rhoE_(ic , j) = p_(ic , j) / (gamma_ - 1);
// 		}
		
// 	}
// }


void Problem::correctInlet()
{
	if (DEBUG)
	{
		std::cout << "Correcting inlet boundary conditions " << std::endl;
	}

	double Vinf = Minf_ * cInf_;
	double Riem1_inf = Vinf + 2/(gamma_-1)*cInf_;

	for (unsigned int i = 0; i < mesh_.Nc_; ++i)
	{
		unsigned int ic = i + 2;
		
		Riem1_(ic , 2) = Riem1_inf;
		Riem2_(ic , 3) = Umag_(ic , 3) - 2/(gamma_-1)*c_(ic , 3);
		Riem2_(ic , 2) = Riem2_(ic , 3);

		Umag_(ic , 2) = 0.5 * (Riem1_(ic , 2) + Riem2_(ic , 2));
		U_(ic , 2) = Umag_(ic , 2);
		c_(ic , 2) = 0.25 * (gamma_ - 1)*(Riem1_(ic, 2) - Riem2_(ic , 2));
		M_(ic , 2) = Umag_(ic , 2)/c_(ic, 2);
		p_(ic , 2) = pInf_/pow((1 + 0.5*(gamma_ - 1) * M_(ic,2) * M_(ic,2)), gamma_/(gamma_-1));
		rho_(ic , 2) = gamma_ * p_(ic , 2)/c_(ic , 2)/c_(ic , 2);
	}

	if (DEBUG)
	{
		std::cout << "U velocity " << std::endl;
		U_.print();

		std::cout << "Density  " << std::endl;
		rho_.print();

		std::cout << "Pressure  " << std::endl;
		p_.print();

		std::cout << "Speed of sound  " << std::endl;
		c_.print();	
	}
}

void Problem::correctWall()
{
	if (DEBUG)
	{
		std::cout << "Correcting wall boundary conditions " << std::endl;
	}

	for (unsigned int j = 0; j < mesh_.Mc_; ++j)
	{
		unsigned int jc = j + 1;
		unsigned int& Icmax = mesh_.Icmax_;

		rho_(Icmax - 1 , jc)  =   rho_(Icmax - 2, jc);	
		U_(Icmax - 1 , jc)    =   U_(Icmax - 2, jc);	
		V_(Icmax - 1 , jc)    = - V_(Icmax - 2, jc);	
		p_(Icmax - 1 , jc)    =   p_(Icmax - 2, jc);
		rhoE_(Icmax - 1 , jc) =   rhoE_(Icmax - 2, jc);	

		rhoU_(Icmax - 1 , jc) = rho_(Icmax - 1 , jc) * U_(Icmax - 1 , jc);	
		rhoV_(Icmax - 1 , jc) = rho_(Icmax - 1 , jc) * V_(Icmax - 1 , jc);	

		rho_(1 , jc)  =   rho_(2 , jc);	
		rhoE_(1 , jc) =   rhoE_(2 , jc);	
		U_(1 , jc)    =   U_(2 , jc);	
		V_(1 , jc)    = - V_(2 , jc);	
		p_(1 , jc)    =   p_(2 , jc);
		rhoU_(1 , jc) = rho_(1 , jc) * U_(1 , jc);	
		rhoV_(1 , jc) = rho_(1 , jc) * V_(1 , jc);	
	}
}

// void Problem::correctPeriodic()
// {
	
// }

void Problem::correctOutlet()
{
	for (unsigned int i = 0; i < Nci_; ++i)
	{
		unsigned int ic = i + 2;
		unsigned int Jmax = mesh_.Jcmax_ - 2; // Last cell

		p_(ic , Jmax ) = pInf_;   // if subsonic!!!!  

		rho_(ic , Jmax ) = rho_(ic , Jmax - 1);

		Riem1_(ic , Jmax ) = Riem1_(ic , Jmax - 1); 

		//V_(ic , Jmax - 1) = V_(ic , Jmax ); 

		rhoU_(ic , Jmax ) = rhoU_(ic , Jmax - 1);

		rhoV_(ic , Jmax ) = rhoV_(ic , Jmax - 1);

		rhoE_(ic , Jmax ) = rhoE_(ic , Jmax - 1);
	}
}

void Problem::correctProperties()
{
	for (int i = 0; i < Nci_; ++i)
	{
		for (int j = 0; j < Mci_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;
			
			Umag_(ic , jc) = sqrt(U_(ic , jc) + V_(ic , jc));
			p_(ic , jc) = (gamma_ - 1) * (rhoE_(ic , jc) - 0.5 * rho_(ic , jc) * Umag_(ic , jc) * Umag_(ic , jc) );
			c_(ic , jc) = sqrt( gamma_ * p_(ic , jc) * rho_(ic , jc));
		}
	}
}

// void Problem::correctBoundaryConditions()
// {
// 	if (DEBUG)
// 	{
// 		std::cout << "Correcting boundary conditions " << std::endl;
// 	}

// 	correctInlet();
// 	correctOutlet();
// 	correctWall();
// }

// Matrix<double> Problem::D(const Matrix<double>& variable)
// {
// 	if (DEBUG)
// 	{
// 		std::cout << "Calculating dissipation " << endl;
// 	}
// 	const unsigned int Nc = mesh_.Nc_;
// 	const unsigned int Mc = mesh_.Mc_;


// 	const double nu2 = 0.5;

// 	// The face area must be corrected 
// 	Matrix<double> yFacesTop = mesh_.yFacesTop_;
// 	Matrix<double> yFacesBottom = mesh_.yFacesBottom_;

// 	Matrix<double> xFacesRight = mesh_.xFacesRight_;
// 	Matrix<double> xFacesLeft = mesh_.xFacesLeft_;

// 	Matrix<double> sCsi2(Nc + 4, Mc + 4);
// 	Matrix<double> sEta2(Nc + 4, Mc + 4);
// 	Matrix<double> lambdaCsi(Nc + 4 , Mc + 4);
// 	Matrix<double> lambdaEta(Nc + 4 , Mc + 4);

// 	Matrix<double> Dvariable(Nc + 4 , Mc + 4);

// 	for (unsigned int i = 0; i < Nc; ++i)
// 	{
// 		std::cout << "i " << i << std::endl;
// 		for (unsigned int j = 0; j < Mc; ++j)
// 		{
// 			unsigned int ic = i + 2;
// 			unsigned int jc = j + 2;

// 			std::cout << "j " << j << std::endl;

// 			double sCsi2_ij_num , sCsi2_ij_den;
// 			double sEta2_ij_num , sEta2_ij_den;

// 			// delta2csi already adjust the indexes to access the fields
// 			sCsi2_ij_num = delta2Csi(p_ , i , j);
// 			sCsi2_ij_den = p_(ic , jc + 1) + p_(ic , jc) + p_(ic , jc - 1);
// 			double sCsi2_ij = sCsi2_ij_num/sCsi2_ij_den;


// 			sEta2_ij_num = delta2Eta(p_ , i , j);
// 			sEta2_ij_den = p_(ic , jc + 1) + p_(ic , jc) + p_(ic , jc - 1);
// 			double sEta2_ij = sEta2_ij_num/sEta2_ij_den;


// 			sCsi2(ic , jc) = sCsi2_ij; 
// 			sEta2(ic , jc) = sEta2_ij; 

// 			lambdaCsi(ic , jc) = std::abs(V_(ic , jc)) + c_(ic , jc);
// 			lambdaEta(ic , jc) = std::abs(U_(ic , jc)) + c_(ic , jc);
// 		}
// 	}



// 	for (unsigned int i = 0; i < Nc; ++i)
// 	{
// 		for (unsigned int j = 0; j < Mc; ++j)
// 		{
// 			unsigned int ic = i + 2;
// 			unsigned int jc = j + 2;

// 			double sCsi2Top = interpolateTop(sCsi2 , i , j);
// 			double sCsi2Bottom = interpolateBottom(sCsi2 , i , j);
// 			double sCsi2Right = interpolateRight(sCsi2 , i , j);
// 			double sCsi2Left = interpolateLeft(sCsi2 , i , j);

// 			double sEta2Top = interpolateTop(sEta2 , i , j);
// 			double sEta2Bottom = interpolateBottom(sEta2 , i , j);
// 			double sEta2Right = interpolateRight(sEta2 , i , j);
// 			double sEta2Left = interpolateLeft(sEta2 , i , j);

// 			double lambdaCsiTop = interpolateTop(lambdaCsi , i , j);
// 			double lambdaCsiBottom = interpolateBottom(lambdaCsi , i , j);
// 			double lambdaCsiRight = interpolateRight(lambdaCsi , i , j);
// 			double lambdaCsiLeft = interpolateLeft(lambdaCsi , i , j);

// 			double lambdaEtaTop = interpolateTop(lambdaEta , i , j);
// 			double lambdaEtaBottom = interpolateBottom(lambdaEta , i , j);
// 			double lambdaEtaRight = interpolateRight(lambdaEta , i , j);
// 			double lambdaEtaLeft = interpolateLeft(lambdaEta , i , j);

// 			double deltaCsi2Var = sCsi2Top * yFacesTop(i , j) * lambdaCsiTop * (variable(ic + 1 , jc) - variable(ic , jc)) 
// 								 - sCsi2Bottom * yFacesBottom(i , j) * lambdaCsiBottom * (variable(ic , jc) - variable(ic - 1 , jc));
	
// 			double deltaEta2Var =  sEta2Right * xFacesRight(i , j) * lambdaEtaRight * (variable(ic , jc + 1) - variable(ic , jc)) 
// 							 	 - sEta2Left * xFacesLeft(i , j) * lambdaEtaLeft * (variable(ic , jc) - variable(ic , jc - 1));

// 			double deltaCsi4Var = sCsi2Top * yFacesTop(i , j) * lambdaCsiTop * ( delta2Csi(variable , i + 1, j) - delta2Csi(variable , i , j) ) 
// 								 - sCsi2Bottom * yFacesBottom(i , j) * lambdaCsiBottom * ( delta2Eta(variable , i , j) - delta2Eta(variable , i - 1, j) );

// 			double deltaEta4Var = sEta2Top * yFacesTop(i , j) * lambdaEtaTop * ( delta2Eta(variable , i , j + 1) - delta2Eta(variable , i , j) ) 
// 								 - sEta2Bottom * yFacesBottom(i , j) * lambdaEtaBottom * ( delta2Eta(variable , i , j ) - delta2Eta(variable , i , j - 1) );

			
// 			Dvariable(i , j)  = (deltaCsi2Var + deltaEta2Var) - (deltaCsi4Var + deltaEta4Var);
			
// 		}
// 	}

// 	return Dvariable;
// }

// void Problem::D(const Matrix<double>& flux_f , const Matrix<double>& flux_g)
// {
// 	if (DEBUG)
// 	{
// 		std::cout << "Calculating dissipation " << endl;
// 	}
// 	const unsigned int Nc = mesh_.Nc_;
// 	const unsigned int Mc = mesh_.Mc_;


// 	const double nu2 = 0.5;

// 	// The face area must be corrected 
// 	Matrix<double> yFacesTop = mesh_.yFacesTop_;
// 	Matrix<double> yFacesBottom = mesh_.yFacesBottom_;

// 	Matrix<double> xFacesRight = mesh_.xFacesRight_;
// 	Matrix<double> xFacesLeft = mesh_.xFacesLeft_;

// 	Matrix<double> sCsi2(Nc + 4, Mc + 4);
// 	Matrix<double> sEta2(Nc + 4, Mc + 4);
// 	Matrix<double> lambdaCsi(Nc + 4 , Mc + 4);
// 	Matrix<double> lambdaEta(Nc + 4 , Mc + 4);


// 	for (unsigned int i = 0; i < Nc; ++i)
// 	{
// 		std::cout << "i " << i << std::endl;
// 		for (unsigned int j = 0; j < Mc; ++j)
// 		{
// 			unsigned int ic = i + 2;
// 			unsigned int jc = j + 2;

// 			std::cout << "j " << j << std::endl;

// 			double sCsi2_ij_num , sCsi2_ij_den;
// 			double sEta2_ij_num , sEta2_ij_den;

// 			// delta2csi already adjust the indexes to access the fields
// 			sCsi2_ij_num = delta2Csi(p_ , i , j);
// 			sCsi2_ij_den = p_(ic , jc + 1) + p_(ic , jc) + p_(ic , jc - 1);
// 			double sCsi2_ij = sCsi2_ij_num/sCsi2_ij_den;


// 			sEta2_ij_num = delta2Eta(p_ , i , j);
// 			sEta2_ij_den = p_(ic , jc + 1) + p_(ic , jc) + p_(ic , jc - 1);
// 			double sEta2_ij = sEta2_ij_num/sEta2_ij_den;


// 			sCsi2(ic , jc) = sCsi2_ij; 
// 			sEta2(ic , jc) = sEta2_ij; 

// 			lambdaCsi(ic , jc) = std::abs(V_(ic , jc)) + c_(ic , jc);
// 			lambdaEta(ic , jc) = std::abs(U_(ic , jc)) + c_(ic , jc);
// 		}
// 	}



// 	for (unsigned int i = 0; i < Nc; ++i)
// 	{
// 		for (unsigned int j = 0; j < Mc; ++j)
// 		{
// 			unsigned int ic = i + 2;
// 			unsigned int jc = j + 2;

// 			double sCsi2Top = interpolateTop(sCsi2 , i , j);
// 			double sCsi2Bottom = interpolateBottom(sCsi2 , i , j);
// 			double sCsi2Right = interpolateRight(sCsi2 , i , j);
// 			double sCsi2Left = interpolateLeft(sCsi2 , i , j);

// 			double sEta2Top = interpolateTop(sEta2 , i , j);
// 			double sEta2Bottom = interpolateBottom(sEta2 , i , j);
// 			double sEta2Right = interpolateRight(sEta2 , i , j);
// 			double sEta2Left = interpolateLeft(sEta2 , i , j);

// 			double lambdaCsiTop = interpolateTop(lambdaCsi , i , j);
// 			double lambdaCsiBottom = interpolateBottom(lambdaCsi , i , j);
// 			double lambdaCsiRight = interpolateRight(lambdaCsi , i , j);
// 			double lambdaCsiLeft = interpolateLeft(lambdaCsi , i , j);

// 			double lambdaEtaTop = interpolateTop(lambdaEta , i , j);
// 			double lambdaEtaBottom = interpolateBottom(lambdaEta , i , j);
// 			double lambdaEtaRight = interpolateRight(lambdaEta , i , j);
// 			double lambdaEtaLeft = interpolateLeft(lambdaEta , i , j);

// 			double deltaCsi2Rho = sCsi2Top * yFacesTop(i , j) * lambdaCsiTop * (rho_(ic + 1 , jc) - rho_(ic , jc)) 
// 								 - sCsi2Bottom * yFacesBottom(i , j) * lambdaCsiBottom * (rho_(ic , jc) - rho_(ic - 1 , jc));
			
// 			double deltaCsi2RhoU = sCsi2Top * yFacesTop(i , j) * lambdaCsiTop * (rhoU_(ic + 1 , jc) - rhoU_(ic , jc)) 
// 								 - sCsi2Bottom * yFacesBottom(i , j) * lambdaCsiBottom * (rho_(ic , jc) - rho_(ic - 1 , jc));
			
// 			double deltaCsi2RhoV = sCsi2Top * yFacesTop(i , j) * lambdaCsiTop * (rhoV_(ic + 1 , jc) - rhoV_(ic , jc)) 
// 							     - sCsi2Bottom * yFacesBottom(i , j) * lambdaCsiBottom * (rho_(ic , jc) - rho_(ic - 1 , jc));
			
// 			double deltaCsi2RhoE = sCsi2Top * yFacesTop(i , j) * lambdaCsiTop * (rhoE_(ic + 1 , jc) - rhoE_(ic , jc)) 
// 								 - sCsi2Bottom * yFacesBottom(i , j) * lambdaCsiBottom * (rho_(ic , jc) - rho_(ic - 1 , jc));




// 			double deltaEta2Rho =  sEta2Right * xFacesRight(i , j) * lambdaEtaRight * (rho_(ic , jc + 1) - rho_(ic , jc)) 
// 							 	 - sEta2Left * xFacesLeft(i , j) * lambdaEtaLeft * (rho_(ic , jc) - rho_(ic , jc - 1));

// 			double deltaEta2RhoU = sEta2Right * xFacesRight(i , j) * lambdaEtaRight * (rhoU_(ic , jc + 1) - rhoU_(ic , jc)) 
// 								 - sEta2Left * xFacesLeft(i , j) * lambdaEtaLeft * (rho_(ic , jc) - rho_(ic , jc - 1));

// 			double deltaEta2RhoV = sEta2Right * xFacesRight(i , j) * lambdaEtaRight * (rhoV_(ic , jc + 1) - rhoV_(ic , jc)) 
// 								 - sEta2Left * xFacesLeft(i , j) * lambdaEtaLeft * (rho_(ic , jc) - rho_(ic , jc - 1));

// 			double deltaEta2RhoE = sEta2Right * xFacesRight(i , j) * lambdaEtaRight * (rhoE_(ic , jc + 1) - rhoE_(ic , jc)) 
// 								 - sEta2Left * xFacesLeft(i , j) * lambdaEtaLeft * (rho_(ic , jc) - rho_(ic , jc - 1));
			
			


// 			double deltaCsi4Rho = sCsi2Top * yFacesTop(i , j) * lambdaCsiTop * ( delta2Csi(rho_ , i + 1, j) - delta2Csi(rho_ , i , j) ) 
// 								 - sCsi2Bottom * yFacesBottom(i , j) * lambdaCsiBottom * ( delta2Eta(rho_ , i , j) - delta2Eta(rho_ , i - 1, j) );

// 			double deltaCsi4RhoU = sCsi2Top * yFacesTop(i , j) * lambdaCsiTop * ( delta2Csi(rhoU_ , i + 1, j) - delta2Csi(rhoU_ , i , j) ) 
// 								 - sCsi2Bottom * yFacesBottom(i , j) * lambdaCsiBottom * ( delta2Eta(rhoU_ , i , j) - delta2Eta(rhoU_ , i - 1 , j) );

// 			double deltaCsi4RhoV = sCsi2Top * yFacesTop(i , j) * lambdaCsiTop * ( delta2Csi(rhoV_ , i + 1, j) - delta2Csi(rhoV_ , i , j) ) 
// 								 - sCsi2Bottom * yFacesBottom(i , j) * lambdaCsiBottom * ( delta2Eta(rhoV_ , i , j) - delta2Eta(rhoV_ , i - 1, j) );

// 			double deltaCsi4RhoE = sCsi2Top * yFacesTop(i , j) * lambdaCsiTop * ( delta2Csi(rhoE_ , i + 1, j) - delta2Csi(rhoE_ , i , j) ) 
// 								 - sCsi2Bottom * yFacesBottom(i , j) * lambdaCsiBottom * ( delta2Eta(rhoE_ , i , j) - delta2Eta(rhoE_ , i - 1, j) );


// 			double deltaEta4Rho = sEta2Top * yFacesTop(i , j) * lambdaEtaTop * ( delta2Eta(rho_ , i , j + 1) - delta2Eta(rho_ , i , j) ) 
// 								 - sEta2Bottom * yFacesBottom(i , j) * lambdaEtaBottom * ( delta2Eta(rho_ , i , j ) - delta2Eta(rho_ , i , j - 1) );

// 			double deltaEta4RhoU = sEta2Top * yFacesTop(i , j) * lambdaEtaTop * ( delta2Eta(rhoU_ , i + 1, j) - delta2Eta(rhoU_ , i , j) ) 
// 								 - sEta2Bottom * yFacesBottom(i , j) * lambdaEtaBottom * ( delta2Eta(rhoU_ , i , j ) - delta2Eta(rhoU_ , i , j - 1) );

// 			double deltaEta4RhoV = sEta2Top * yFacesTop(i , j) * lambdaEtaTop * ( delta2Eta(rhoV_ , i + 1, j) - delta2Eta(rhoV_ , i , j) ) 
// 								 - sEta2Bottom * yFacesBottom(i , j) * lambdaEtaBottom * ( delta2Eta(rhoV_ , i , j ) - delta2Eta(rhoV_ , i , j - 1) );

// 			double deltaEta4RhoE = sEta2Top * yFacesTop(i , j) * lambdaEtaTop * ( delta2Eta(rhoE_ , i + 1, j) - delta2Eta(rhoE_ , i , j) ) 
// 								 - sEta2Bottom * yFacesBottom(i , j) * lambdaEtaBottom * ( delta2Eta(rhoE_ , i , j ) - delta2Eta(rhoE_ , i , j - 1) );
			

// 			Drho_(i , j)  = (deltaCsi2Rho + deltaEta2Rho) - (deltaCsi4Rho + deltaEta4Rho);
// 			DrhoU_(i , j) = (deltaCsi2RhoU + deltaEta2RhoU) - (deltaCsi4RhoU + deltaEta4RhoU);
// 			DrhoV_(i , j) = (deltaCsi2RhoV + deltaEta2RhoV) - (deltaCsi4RhoV + deltaEta4RhoV);
// 			DrhoE_(i , j) = (deltaCsi2RhoE + deltaEta2RhoE) - (deltaCsi4RhoE + deltaEta4RhoE);
// 		}
// 	}
// }


// Matrix<double> Problem::R(const Matrix<double>& flux_f , const Matrix<double>& flux_g)
// {
// 	const unsigned int Nc = mesh_.Nc_;
// 	const unsigned int Mc = mesh_.Mc_;

// 	Matrix<double>& xFacesLeft = mesh_.xFacesLeft_;
// 	Matrix<double>& xFacesRight = mesh_.xFacesRight_;
// 	Matrix<double>& xFacesTop = mesh_.xFacesTop_;
// 	Matrix<double>& xFacesBottom = mesh_.xFacesBottom_;


// 	Matrix<double>& yFacesLeft = mesh_.yFacesLeft_;
// 	Matrix<double>& yFacesRight = mesh_.yFacesRight_;
// 	Matrix<double>& yFacesTop = mesh_.yFacesTop_;
// 	Matrix<double>& yFacesBottom = mesh_.yFacesBottom_;

// 	// Nc is the number of internal cells
// 	Matrix<double> Rvariable(Nc + 4 , Mc + 4);


// 	for (unsigned int i = 0; i < Nc; ++i)
// 	{
// 		for (unsigned int j = 0; j < Mc; ++j)
// 		{
// 			unsigned int ic = i + 2;
// 			unsigned int jc = j + 2;
			
// 			// Declare fluxes cell by cell at the left,right,top and bottom locations
// 			double flux_f_Left, flux_f_Right, flux_f_Top, flux_f_Bottom;
// 			double flux_g_Left, flux_g_Right, flux_g_Top, flux_g_Bottom;

// 			// Interpolate fields and compute fluxes

// 			// // f fluxes: for rectilinear grids it concerns only left and right. 
// 			// // More general grids require xFacesTop and xFacesTop
// 			// // Example rhoPhi: rhoPhi_f = rhoPhiLeft + rhoPhiRight + rhoPhiBottom + rhoPhiTop
// 			// // The difference between rhoPhi for f and g is the velocity component

// 			// f
// 			flux_f_Left       = interpolateLeft(flux_f  , i , j) * xFacesLeft(i , j);
// 			flux_f_Right      = interpolateRight(flux_f , i , j) * xFacesRight(i , j);
// 			flux_f_Top        = interpolateTop(flux_f , i , j) * xFacesTop(i , j);
// 			flux_f_Bottom     = interpolateBottom(flux_f , i , j) * xFacesBottom(i , j);

// 			std::cout << "i " << i << std::endl;
// 			std::cout << "j " << j << std::endl;
// 			std::cout << "flux_f_Left " << flux_f_Left << std::endl;
// 			std::cout << "flux_f_Right " << flux_f_Right << std::endl;
// 			std::cout << "flux_f_Top " << flux_f_Top << std::endl;
// 			std::cout << "flux_f_Bottom " << flux_f_Bottom << std::endl;

// 			// g
// 			flux_g_Left       = interpolateLeft(flux_g  , i , j) * yFacesLeft(i , j);
// 			flux_g_Right      = interpolateRight(flux_g , i , j) * yFacesRight(i , j);
// 			flux_g_Top        = interpolateTop(flux_g , i , j) * yFacesTop(i , j);
// 			flux_g_Bottom     = interpolateBottom(flux_g , i , j) * yFacesBottom(i , j);
					

// 			// Compute R-component of fluxes
// 			Rvariable(ic,jc) = (flux_f_Right - flux_f_Left + flux_f_Top - flux_f_Bottom) -
// 							   (flux_g_Right - flux_g_Left + flux_g_Top - flux_g_Bottom);
// 		}
// 	}

// 	return Rvariable;


// }

// void Problem::R()
// {
	// if (DEBUG)
	// {
	// 	std::cout << "Calculating residuals " << endl;
	// }
	// const unsigned int Nc = mesh_.Nc_;
	// const unsigned int Mc = mesh_.Mc_;

	// Matrix<double>& xFacesLeft = mesh_.xFacesLeft_;
	// Matrix<double>& xFacesRight = mesh_.xFacesRight_;
	// Matrix<double>& xFacesTop = mesh_.xFacesTop_;
	// Matrix<double>& xFacesBottom = mesh_.xFacesBottom_;


	// Matrix<double>& yFacesLeft = mesh_.yFacesLeft_;
	// Matrix<double>& yFacesRight = mesh_.yFacesRight_;
	// Matrix<double>& yFacesTop = mesh_.yFacesTop_;
	// Matrix<double>& yFacesBottom = mesh_.yFacesBottom_;


	// for (unsigned int i = 0; i < Nc; ++i)
	// {
	// 	for (unsigned int j = 0; j < Mc; ++j)
	// 	{
	// 		unsigned int ic = i + 2;
	// 		unsigned int jc = j + 2;
			
	// 		 U_(ic,jc) = rhoU_(ic,jc)/rho_(ic,jc);
	// 		 V_(ic,jc) = rhoV_(ic,jc)/rho_(ic,jc);
	// 		 rhoU2_(ic,jc) = rhoU_(ic,jc)*U_(ic,jc);
	// 		 rhoV2_(ic,jc) = rhoV_(ic,jc)*V_(ic,jc);
	// 		 rhoUV_(ic,jc) = rhoU_(ic,jc)*V_(ic,jc);
	// 		 rhoHu_(ic,jc) = rhoE_(ic,jc)*U_(ic,jc) + p_(ic,jc);
	// 		 rhoHv_(ic,jc) = rhoE_(ic,jc)*V_(ic,jc) + p_(ic,jc);

	// 		// Declare fluxes cell by cell at the left,right,top and bottom locations
	// 		double rhoUPhiLeft, rhoUPhiRight, rhoUPhiTop, rhoUPhiBottom;

	// 		double rhoVPhiLeft, rhoVPhiRight, rhoVPhiTop, rhoVPhiBottom;

	// 		double rhoU2PhiLeft, rhoU2PhiRight, rhoU2PhiTop, rhoU2PhiBottom;

	// 		double pLeft_f, pRight_f, pTop_f, pBottom_f;
			
	// 		double rhoUVPhiLeft_f, rhoUVPhiRight_f, rhoUVPhiTop_f, rhoUVPhiBottom_f;

			
	// 		double rhoV2PhiLeft, rhoV2PhiRight, rhoV2PhiTop, rhoV2PhiBottom;
			
	// 		double pLeft_g, pRight_g, pTop_g, pBottom_g;

	// 		double rhoUVPhiLeft_g, rhoUVPhiRight_g, rhoUVPhiTop_g, rhoUVPhiBottom_g;

	// 		double rhoHuPhiLeft, rhoHuPhiRight, rhoHuPhiTop, rhoHuPhiBottom;

	// 		double rhoHvPhiLeft, rhoHvPhiRight, rhoHvPhiTop, rhoHvPhiBottom;


	// 		// Interpolate fields and compute fluxes

	// 		// // f fluxes: for rectilinear grids it concerns only left and right. 
	// 		// // More general grids require xFacesTop and xFacesTop
	// 		// // Example rhoPhi: rhoPhi_f = rhoPhiLeft + rhoPhiRight + rhoPhiBottom + rhoPhiTop
	// 		// // The difference between rhoPhi for f and g is the velocity component

	// 		// f
	// 		rhoUPhiLeft       = interpolateLeft(rhoU_  , i , j) * xFacesLeft(i , j);
	// 		rhoUPhiRight      = interpolateRight(rhoU_ , i , j) * xFacesRight(i , j);
	// 		rhoUPhiTop        = interpolateTop(rhoU_ , i , j) * xFacesTop(i , j);
	// 		rhoUPhiBottom     = interpolateBottom(rhoU_ , i , j) * xFacesBottom(i , j);

	// 		std::cout << "i " << i << std::endl;
	// 		std::cout << "j " << j << std::endl;
	// 		std::cout << "rhoUPhiLeft " << rhoUPhiLeft << std::endl;
	// 		std::cout << "rhoUPhiRight " << rhoUPhiRight << std::endl;
	// 		std::cout << "rhoUPhiTop " << rhoUPhiTop << std::endl;
	// 		std::cout << "rhoUPhiBottom " << rhoUPhiBottom << std::endl;

			

	// 		rhoU2PhiLeft      = interpolateLeft(rhoU2_   , i , j) * xFacesLeft(i , j);
	// 		rhoU2PhiRight     = interpolateRight(rhoU2_  , i , j) * xFacesRight(i , j);
	// 		rhoU2PhiTop       = interpolateTop(rhoU2_    , i , j) * xFacesTop(i , j);
	// 		rhoU2PhiBottom    = interpolateBottom(rhoU2_ , i , j) * xFacesBottom(i , j);
			
	// 		rhoUVPhiLeft_f      = interpolateLeft(rhoUV_   , i , j) * xFacesLeft(i , j); 
	// 		rhoUVPhiRight_f     = interpolateRight(rhoUV_  , i , j) * xFacesRight(i , j);
	// 		rhoUVPhiTop_f       = interpolateTop(rhoUV_    , i , j) * xFacesTop(i , j);
	// 		rhoUVPhiBottom_f    = interpolateBottom(rhoUV_ , i , j) * xFacesBottom(i , j);
			
	// 		rhoHuPhiLeft      = interpolateLeft(rhoHu_  , i , j) * xFacesLeft(i , j);
	// 		rhoHuPhiRight     = interpolateRight(rhoHu_ , i , j) * xFacesRight(i , j);
	// 		rhoHuPhiTop       = interpolateTop(rhoHu_   , i , j) * xFacesTop(i , j);
	// 		rhoHuPhiBottom    = interpolateBottom(rhoHu_  , i , j) * xFacesBottom(i , j);
			
	// 		pLeft_f          = interpolateLeft(p_  , i , j) * xFacesLeft(i , j);
	// 		pRight_f         = interpolateRight(p_   , i , j) * xFacesRight(i , j);
	// 		pTop_f           = interpolateTop(p_  , i , j) * xFacesTop(i , j);
	// 		pBottom_f        = interpolateBottom(p_  , i , j) * xFacesBottom(i , j);
			
			
	// 		// g
	// 		rhoVPhiLeft      = interpolateLeft  (rhoV_    , i , j) * yFacesLeft(i , j);
	// 		rhoVPhiRight     = interpolateRight (rhoV_   , i , j) * yFacesRight(i , j);
	// 		rhoVPhiTop       = interpolateTop   (rhoV_     , i , j) * yFacesTop(i , j);
	// 		rhoVPhiBottom    = interpolateBottom(rhoV_  , i , j) * yFacesBottom(i , j);
	// 		std::cout << "rhoVPhiLeft " << rhoVPhiLeft << std::endl;
	// 		std::cout << "rhoVPhiRight " << rhoVPhiRight << std::endl;
	// 		std::cout << "rhoVPhiTop " << rhoVPhiTop << std::endl;
	// 		std::cout << "rhoVPhiBottom " << rhoVPhiBottom << std::endl;
			
	// 		rhoV2PhiLeft      = interpolateLeft  (rhoV2_    , i , j) * yFacesLeft(i , j);
	// 		rhoV2PhiRight     = interpolateRight (rhoV2_   , i , j) * yFacesRight(i , j);
	// 		rhoV2PhiTop       = interpolateTop   (rhoV2_     , i , j) * yFacesTop(i , j);
	// 		rhoV2PhiBottom    = interpolateBottom(rhoV2_  , i , j) * yFacesBottom(i , j);
			
	// 		rhoUVPhiLeft_g      = interpolateLeft  (rhoUV_    , i , j) * yFacesLeft(i , j);
	// 		rhoUVPhiRight_g     = interpolateRight (rhoUV_   , i , j) * yFacesRight(i , j);
	// 		rhoUVPhiTop_g       = interpolateTop   (rhoUV_     , i , j) * yFacesTop(i , j);
	// 		rhoUVPhiBottom_g    = interpolateBottom(rhoUV_  , i , j) * yFacesBottom(i , j);
			
	// 		pLeft_g      = interpolateLeft  (p_    , i , j) * yFacesLeft(i , j);
	// 		pRight_g     = interpolateRight (p_	   , i , j) * yFacesRight(i , j);
	// 		pTop_g       = interpolateTop   (p_     , i , j) * yFacesTop(i , j);
	// 		pBottom_g    = interpolateBottom(p_  , i , j) * yFacesBottom(i , j);
			
	// 		rhoHvPhiLeft     = interpolateLeft  (p_    , i , j) * yFacesLeft(i , j);
	// 		rhoHvPhiRight    = interpolateRight  (p_    , i , j) * yFacesRight(i , j);
	// 		rhoHvPhiTop      = interpolateTop  (p_   , i , j) * yFacesTop(i , j);
	// 		rhoHvPhiBottom   = interpolateBottom  (p_    , i , j) * yFacesBottom(i , j);
			

	// 		// Compute R-component of residuals
	// 		// R(i,j) = sum_l (f_l * deltaY_l - g_l * deltaX_l)
	// 		Rrho_(ic,jc) = (rhoUPhiRight - rhoUPhiLeft + rhoUPhiTop - rhoUPhiBottom) - (rhoVPhiRight - rhoVPhiLeft + rhoVPhiTop - rhoVPhiBottom);
	// 		RrhoU_(ic,jc) = (rhoU2PhiRight + pRight_f) - (rhoU2PhiLeft + pLeft_f) + (rhoU2PhiTop + pTop_f) - (rhoU2PhiBottom + pBottom_f) -
	// 		   				(rhoUVPhiRight_g - rhoUVPhiLeft_g + rhoUVPhiTop_g - rhoUVPhiBottom_g);
	// 		RrhoV_(ic,jc) = (rhoV2PhiRight + pRight_g) - (rhoV2PhiLeft + pLeft_g) + (rhoV2PhiTop + pTop_g) - (rhoV2PhiBottom + pBottom_g) -
	// 		   				(rhoUVPhiRight_f - rhoUVPhiLeft_f + rhoUVPhiTop_f - rhoUVPhiBottom_f);;
	// 		RrhoE_(ic,jc) = (rhoHuPhiRight - rhoHuPhiLeft + rhoHuPhiTop - rhoHuPhiBottom) - (rhoHvPhiRight - rhoHvPhiLeft + rhoHvPhiTop - rhoHvPhiBottom);
	// 	}
	// }
//}


// void Problem::RungeKutta(Matrix<double>& variable, Matrix<double>& R , Matrix<double> D0 , int i , int j )
// {
// 	double alpha1 = 0.25, alpha2 = 0.5 , alpha3 = 1/3, alpha4 = 1;
// 	double dt = 1e-5;

// 	Matrix<double>& area = mesh_.area_;

// 	int ic = i + 2;
// 	int jc = j + 2;

// 	double Aij = area(ic,jc);

// 	// variable(ic , jc) = variable(ic , jc) - alpha1 * dt/Aij * (R(ic , jc) - D(ic , jc));

// }

// void Problem::solve()
// {
// 	// const unsigned int Nc = mesh_.Nc_ + 1;
// 	// const unsigned int Mc = mesh_.Mc_ + 1;

// 	// Matrix<double> Rrho(Nc + 3 , Mc + 3);
// 	// Matrix<double> RrhoU(Nc + 3 , Mc + 3);
// 	// Matrix<double> RrhoV(Nc + 3 , Mc + 3);
// 	// Matrix<double> RrhoE(Nc + 3 , Mc + 3);

// 	// Rrho = R(rhoU_ , rhoV_)
// 	// RrhoU = R(rhoU_ , rhoV_)

// 	// Matrix<double> Frho(Nc + 3 , Mc + 3);
// 	// Matrix<double> FrhoU(Nc + 3 , Mc + 3);
// 	// Matrix<double> FrhoV(Nc + 3 , Mc + 3);
// 	// Matrix<double> FrhoE(Nc + 3 , Mc + 3);

	
// 	// int n = 0;
// 	// while(n < 1)
// 	// {	

// 	// 	R();
// 	// 	D();
// 	// 	for (unsigned int i = 0; i < Nc; ++i)
// 	// 	{
// 	// 		for (unsigned int j = 0; j < Mc; ++j)
// 	// 		{
// 	// 			unsigned int ic = i + 2;
// 	// 			unsigned int jc = j + 2;
// 	// 			Frho(ic , jc) = Rrho_(ic , jc) - Drho_(ic , jc);
// 	// 			FrhoU(ic , jc) = RrhoU_(ic , jc) - DrhoU_(ic , jc);
// 	// 			FrhoV(ic , jc) = RrhoV_(ic , jc) - DrhoV_(ic , jc);
// 	// 			FrhoE(ic , jc) = RrhoE_(ic , jc) - DrhoE_(ic , jc);

// 	// 		}
// 	// 	}

// 	// 	if (DEBUG)
// 	// 	{
// 	// 		std::cout << "Mass flux Frho \n " << std::endl;
// 	// 		Frho.print(); 
// 	// 		std::cout << "Momentum flux FrhoU \n " << std::endl;
// 	// 		FrhoU.print(); 
// 	// 		std::cout << "Momentum flux FrhoV \n " << std::endl;
// 	// 		FrhoV.print(); 
// 	// 		std::cout << "Momentum flux FrhoE \n " << std::endl;
// 	// 		FrhoE.print(); 
// 	// 	}

// 	// 	n++;

// 	// }
// }