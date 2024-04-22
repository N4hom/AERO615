#include "Matrix.hpp"
#include "Mesh.hpp"
#include "Variable.hpp"
#include "debug.hpp"
#include <math.h>
#include <algorithm>
#include <ostream>
#include <fstream>

class Problem
{

public:
	Mesh mesh_;
	double dt_;
	double CFL_ = 2.8;
	const double nu2_ = 0;
	const double nu4_ = 0.001;
	unsigned int Nci_ , Mci_;
	unsigned int Nc_ , Mc_;

	Matrix<double> U_, V_ , p_, T_ , Umag_;							   // Primitive variables 
	Variable        rho_, rhoU_,rhoV_, rhoE_;  // State vector components
	Matrix<double> c_;
	Matrix<double> M_;
	Matrix<double> Riem1_;
	Matrix<double> Riem2_;
	Matrix<double> s2Csi_;
	Matrix<double> s2Eta_;
	Matrix<Coefficients> s2_;
	Matrix<Coefficients> s4_;
	Matrix<Coefficients> lambda_;

	double gamma_ = 1.4;
	double pInf_  = 101325;;
	double pRatio = 0.9;
	double Minf_  = 0.3;
	double alphaInf_ = 0;
	double rhoInf_   = 1;
	double cInf_     = sqrt(gamma_*pInf_/rhoInf_);
	double Uinf_  = Minf_ * cInf_;  

	std::ofstream outputFile;

	
	double delta2Csi(const Matrix<double>& matrix, unsigned int i, unsigned int j);
	double delta2Eta(const Matrix<double>& matrix, unsigned int i, unsigned int j);
	
	double interpolateLeft  (Matrix<double>& flux , unsigned int i, unsigned int j) const;
	double interpolateRight (Matrix<double>& flux , unsigned int i, unsigned int j) const;
	double interpolateTop   (Matrix<double>& flux , unsigned int i, unsigned int j) const;
	double interpolateBottom(Matrix<double>& flux , unsigned int i, unsigned int j) const;

	void initialize();
	void computeFluxes();
	void computeSwitches();
	void correctBoundaryConditions();
	void correctInlet();
	void correctWall();
	void correctPeriodic();
	void correctOutlet();
	void correctProperties();
	void correctPropertiesij(unsigned int i , unsigned int j);
	void correctTimeStep();
	void R();
	// Matrix<double> R(const Matrix<double>& flux_f , const Matrix<double>& flux_g);
	// Matrix<double> D(const Matrix<double>& variable);
	//void D(const Matrix<double>& flux_f , const Matrix<double>& flux_g);
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
U_(Nc_ , Mc_),
V_(Nc_ , Mc_),
p_(Nc_ , Mc_),
T_(Nc_ , Mc_),
c_(Nc_ , Mc_),
M_(Nc_ , Mc_),
Umag_(Nc_ , Mc_),
Riem1_(Nc_ , Mc_),
Riem2_(Nc_ , Mc_),
s2_(Nc_ , Mc_),
s4_(Nc_ , Mc_),
s2Csi_(Nc_ , Mc_),
s2Eta_(Nc_ , Mc_),
lambda_(Nc_ , Mc_),
outputFile("Residuals"),
rho_("rho"   , Nc_ , Mc_, mesh_, U_, V_, c_ , p_ , s2_ , s4_),  // N - 1 being the number of internal cells
rhoU_("rhoU" , Nc_ , Mc_, mesh_, U_, V_, c_ , p_ , s2_ , s4_),
rhoV_("rhoV" , Nc_ , Mc_, mesh_, U_, V_, c_ , p_ , s2_ , s4_),
rhoE_("rhoE" , Nc_ , Mc_, mesh_, U_, V_, c_ , p_ , s2_ , s4_)
{

	initialize();
	rho_.print();
	rhoU_.print();
	rhoV_.print();
	rhoE_.print();
	p_.print();
	c_.print();
	U_.print();
	V_.print();

	


	correctInlet();
	correctWall();
	correctOutlet();
	rho_.print();
	rhoU_.print();
	rhoV_.print();
	rhoE_.print();

	computeFluxes();

	if(DEBUG)
	{	
		std::cout<< "rho flux f:" << std::endl;
		rho_.flux_f().print();
		std::cout<< "rho flux g:" << std::endl;
		rho_.flux_g().print();
		std::cout<< "rhoU flux f:" << std::endl;
		rhoU_.flux_f().print();
		std::cout<< "rhoU flux g:" << std::endl;
		rhoU_.flux_g().print();
		std::cout<< "rhoV flux f:" << std::endl;
		rhoV_.flux_f().print();
		std::cout<< "rhoV flux g:" << std::endl;
		rhoV_.flux_g().print();
		std::cout<< "rhoE flux f:" << std::endl;
		rhoE_.flux_f().print();
		std::cout<< "rhoE flux g:" << std::endl;
		rhoE_.flux_g().print();
	}

	rho_.print();
	rhoU_.print();
	rhoV_.print();
	rhoE_.print();
	p_.print();
	c_.print();
	U_.print();
	V_.print();

	// It should be executed inside computeDissipation() once every time step
	computeSwitches();
	rho_.computeResidual();
	rhoU_.computeResidual();
	rhoV_.computeResidual();
	rhoE_.computeResidual();
	rho_.computeDissipation();
	rhoU_.computeDissipation();
	rhoV_.computeDissipation();
	rhoE_.computeDissipation();
	
	
}

Problem::~Problem(){}


double Problem::interpolateLeft(Matrix<double>& flux , unsigned int i, unsigned int j) const
{
	
	return 0.5 * (flux(i     , j - 1) + flux(i , j)) ;
}

double Problem::interpolateRight(Matrix<double>& flux , unsigned int i, unsigned int j) const
{
	return 0.5 * (flux(i     , j + 1) + flux(i , j)) ;
}

double Problem::interpolateTop(Matrix<double>& flux , unsigned int i, unsigned int j)	const
{
	return 0.5 * (flux(i - 1  , j   ) + flux(i , j)) ;
}

double Problem::interpolateBottom(Matrix<double>& flux , unsigned int i, unsigned int j) const
{
	return 0.5 * (flux(i + 1  , j   ) + flux(i , j)) ;
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
	
	return matrix(ic , jc + 1) - 2*matrix(ic , jc) + matrix(ic , jc - 1);
}

void Problem::initialize()
{

	std::cout << "Initializing internal field and ghost cells at outlet" << std::endl;
	// Initialize internal field
	for (unsigned int ic = 1; ic < Nc_ - 1; ++ic)
	{
		for (unsigned int jc = 0; jc < Mc_ ; ++jc)
		{
			// unsigned ic = i + 2;
			// unsigned jc = j + 2;
			rho_(ic , jc) = rhoInf_*pRatio;
			p_(ic , jc) = pInf_*pRatio;
			c_(ic , jc) = cInf_;
			U_(ic , jc) = Minf_ * cInf_;
			rhoU_(ic , jc) = rho_(ic , jc) * U_(ic , jc);
			rhoV_(ic , jc) = rho_(ic , jc) * V_(ic , jc);
			Umag_(ic , jc) = sqrt(U_(ic , jc)*U_(ic , jc) + V_(ic , jc) * V_(ic , jc));
			rhoE_(ic , jc) = p_(ic , jc)/(gamma_ - 1) + 0.5 * rho_(ic , jc) * Umag_(ic , jc)*Umag_(ic , jc);
		}
	}

	// std::cout << "Initializing inlet " << std::endl;
	// //Initialize inlet boundary
	// for (unsigned int ic = 1; ic < Nc_ - 1; ++ic)
	// {
	// 	for (int j = 0 ; j < 3; ++j)
	// 	{
				
	// 		// unsigned ic = i + 2; 
	// 		// I don't use jc here because there's no need to since I'm setting the boundary conditions to
	// 		// the initial values. Only the first two columns of cells are touched.

	// 		rho_(ic , j)  = rhoInf_;
	// 		p_(ic , j)    = pInf_;
	// 		c_(ic , j)    = cInf_;
	// 		U_(ic , j )   = Minf_ * cInf_;
	// 		rhoU_(ic , j) = rho_(ic , j) * U_(ic , j);
	// 		rhoV_(ic , j) = rho_(ic , j) * V_(ic , j);
	// 		Umag_(ic , j) = sqrt(U_(ic , j)*U_(ic , j) + V_(ic , j) * V_(ic , j));
	// 		rhoE_(ic , j) = p_(ic , j) / (gamma_ - 1) + 0.5 * rho_(ic , j) * Umag_(ic , j)*Umag_(ic , j);
	// 	}
		
	// }
	
}

void Problem::computeSwitches()
{
	for (unsigned int i = 0; i < Nci_; ++i)
	{
		for (unsigned int j = 0; j < Mci_; ++j)
		{
			

			unsigned int ic = i + 2;
			unsigned int jc = j + 2;

			s2Csi_(ic , jc) = nu2_ * std::abs(delta2Csi(p_ , i , j))/(p_(ic , jc + 1) + 2*p_(ic , jc) + p_(ic , jc - 1));
			s2Eta_(ic , jc) = nu2_ * std::abs(delta2Eta(p_ , i , j))/(p_(ic , jc + 1) + 2*p_(ic , jc) + p_(ic , jc - 1));

			
		}
			
	}


	// Interpolate switches at cell faces
	// i goes from 1 to Nci_ - 1 because the switches at the north and south boundaries are supposed to be zero
	for (unsigned int i = 0; i < Nci_ ; ++i)
	{
		for (unsigned int j = 0; j < Mci_; ++j)
		{

			unsigned int ic = i + 2;
			unsigned int jc = j + 2;


			s2_(ic , jc).e = interpolateRight(s2Eta_ , ic , jc);
			s2_(ic  , jc).w = interpolateLeft(s2Eta_ , ic , jc);

			if (i == 0)
			{
				s2_(ic , jc).n = 0;
				s2_(ic , jc).s = interpolateBottom(s2Csi_ , ic , jc);
			}
			else if (i == Nci_ - 1)
			{
				s2_(ic , jc).s = 0;
				s2_(ic , jc).n = interpolateTop(s2Csi_ , ic , jc);

			}
			else
			{
				s2_(ic , jc).n = interpolateTop(s2Csi_ , ic , jc);
				s2_(ic , jc).s = interpolateBottom(s2Csi_ , ic , jc);

			}

			
			s4_(ic , jc).e = std::max(0. , nu4_ - s2_(ic , jc).e );
			s4_(ic , jc).w = std::max(0. , nu4_ - s2_(ic , jc).w );

			if (i == 0)
			{
				s4_(ic , jc).n = 0;
				s4_(ic , jc).s = std::max(0. , nu4_ - s2_(ic , jc).s );

			}
			else if (i == Nci_ - 1)
			{

				s4_(ic , jc).n = std::max(0. , nu4_ - s2_(ic , jc).n );
				s4_(ic , jc).s = 0;

			}
			else
			{
				s4_(ic , jc).n = std::max(0. , nu4_ - s2_(ic , jc).n );
				s4_(ic , jc).s = std::max(0. , nu4_ - s2_(ic , jc).s );
			}


			
		}
			
	}


	// Printing
	if (DEBUG)
	{
		// Print s2
		for (unsigned int i = 0; i < Nci_ ; ++i)
		{
			for (unsigned int j = 0; j < Mci_; ++j)
			{

				unsigned int ic = i + 2;
				unsigned int jc = j + 2;

				std::cout << "i " << i << std::endl;
				std::cout << "j " << j << std::endl;
				std::cout << "s2 " << std::endl;
				s2_(ic , jc).print();
				
			}
				
		}

		// Print s4
		for (unsigned int i = 0; i < Nci_ ; ++i)
		{
			for (unsigned int j = 0; j < Mci_; ++j)
			{

				unsigned int ic = i + 2;
				unsigned int jc = j + 2;

				std::cout << "i " << i << std::endl;
				std::cout << "j " << j << std::endl;
				std::cout << "s4 " << std::endl;
				s4_(ic , jc).print();
			}
				
		}
		/* code */
	}

}

void Problem::computeFluxes()
{
	if (DEBUG)
	{
		std::cout << "Computing fluxes " << std::endl;
	}
	
	// Looping over all domain because the fluxes in the ghost cells are needed
	for (unsigned int i = 0; i < Nc_; ++i)
	{
		for (unsigned int j = 0; j < Mc_; ++j)
		{
			// unsigned ic = i + 2;
			// unsigned jc = j + 2;

			rho_.correctFlux_f(i , j);
			rho_.correctFlux_g(i , j);
			
			rhoU_.correctFlux_f(i , j);
			rhoU_.correctFlux_g(i , j);
			
			rhoV_.correctFlux_f(i , j);
			rhoV_.correctFlux_g(i , j);
			
			rhoE_.correctFlux_f(i , j);
			rhoE_.correctFlux_g(i , j);

			// rhoFlux_f(ic , jc)  = rhoU_ (ic , jc);
			// rhoUFlux_f(ic , jc) = rhoU_(ic , jc) * U_(ic , jc) + p_(ic , jc);
			// rhoVFlux_f(ic , jc) = rhoV_(ic , jc) * U_(ic , jc);
			// rhoEFlux_f(ic , jc) = rhoE_(ic , jc) * U_(ic , jc) + p_(ic , jc) * U_(ic , jc);
			
			// rhoFlux_g(ic , jc)  = rhoV_ (ic , jc);
			// rhoUFlux_g(ic , jc) = rhoU_(ic , jc) * V_(ic , jc);
			// rhoVFlux_g(ic , jc) = rhoV_(ic , jc) * V_(ic , jc) + p_(ic , jc);
			// rhoEFlux_g(ic , jc) = rhoE_(ic , jc) * V_(ic , jc) + p_(ic , jc) * V_(ic , jc);
			
		}
	}


	

}


void Problem::solve()
{

	unsigned int N =80;
	unsigned int iter = 0;


	while(iter < N)
	{
		if (DEBUG)
		{
			std::cout << "-------------------------------------------------------------------" << std::endl;
			std::cout << "iter " << iter << std::endl;
			std::cout << "dt " << dt_ << std::endl;
		}
			
		correctTimeStep();
		

		rho_.computeResidual();
		rho_.computeDissipation();
		rhoU_.computeResidual();
		rhoU_.computeDissipation();
		rhoV_.computeResidual();
		rhoV_.computeDissipation();
		rhoE_.computeResidual();
		rhoE_.computeDissipation();
		for (int i = 0; i < Nci_; ++i)
		{
			for (int j = 0; j < Mci_  ; ++j)
			{
				unsigned int ic = i + 2;
				unsigned int jc = j + 2;

				rho_.phiPrev()(ic,jc) = rho_(ic , jc);
				RungeKutta(rho_  , i , j);
				rho_.print();
				RungeKutta(rhoU_  , i , j);
				rhoU_.print();
				RungeKutta(rhoV_  , i , j);
				rhoV_.print();
				RungeKutta(rhoE_  , i , j);
				rhoE_.print();
				std::cout << "\n " << std::endl;
				//double& Aij       = area(ic,jc);
			}
		}

		
		rho_.print();
		rho_.printFlux_f();
		rho_.printFlux_g();
		rhoU_.print();
			// std::cout << "rhoV " << std::endl;
		rhoV_.print();
			// std::cout << "rhoE " << std::endl;
		rhoE_.print();
		std::cout << "U " << std::endl;
		U_.print();
		std::cout << "V " << std::endl;
		V_.print();
		std::cout << "p " << std::endl;
		p_.print();
		std::cout << "c " << std::endl;
		c_.print();
				
		
		correctProperties();
		correctWall();
		correctInlet();
		correctOutlet();

		computeSwitches();
		computeFluxes();

		

		std::cout << "End of time step " << std::endl;
		if (DEBUG)
		{
			outputFile <<  rho_.computeError() << std::endl;

			rho_.print();
			rho_.printFlux_f();
			rho_.printFlux_g();
			rhoU_.print();
			// std::cout << "rhoV " << std::endl;
			rhoV_.print();
			// std::cout << "rhoE " << std::endl;
			rhoE_.print();

			std::cout << "U " << std::endl;
			U_.print();
			std::cout << "V " << std::endl;
			V_.print();
			std::cout << "p " << std::endl;
			p_.print();
			std::cout << "c " << std::endl;
			c_.print();

			// std::cout << "M " << std::endl;
			// M_.print();

		}
		else
		{
			outputFile <<  rho_.computeError() << std::endl;
			rho_.print();
			rhoU_.print();
			rhoV_.print();
			rhoE_.print();

			std::cout << "U " << std::endl;
			U_.print();
			std::cout << "V " << std::endl;
			V_.print();
			std::cout << "p " << std::endl;
			p_.print();
			std::cout << "c " << std::endl;
			c_.print();
		}

		std::cout << "--------------------------------------------------------------------------------" << std::endl;

		iter++;
	}





}


void Problem::RungeKutta(Variable& variable, int i , int j )
{	
	if (DEBUG)
	{std::cout << "\n " << std::endl;
		std::cout << "Runge-Kutta " << "variable " << variable.name() << std::endl;
		std::cout << "i " << i << std::endl;
		std::cout << "j " << j << std::endl;
		// std::cout << "R "  <<  std::endl;
		variable.R().print();
		variable.print() ;
		variable.printFlux_f() ;
		variable.printFlux_g() ;
	}

	const double alpha1 = 0.25, alpha2 = 0.5 , alpha3 = 0.33333, alpha4 = 1.;
	// std::cout << "alpha1 " << alpha1 << std::endl;
	// std::cout << "alpha2 " << alpha2 << std::endl;
	// std::cout << "alpha3 " << alpha3 << std::endl;
	// std::cout << "alpha4 " << alpha4 << std::endl;

	Matrix<double>& area = mesh_.area_;

	int ic = i + 2;
	int jc = j + 2;

	double& Aij       = area(ic,jc);
	double& Rij       = variable.R()(i , j);
	double  Dij0      = variable.D()(i , j);
	double  variable0 = variable(ic , jc);
	double dtByAij = dt_/Aij;
	// At the next RK step the fluxes must be updated to calculate the new residual. It will be changed only the flux and residual at the cell (ic jc) or (i j)
	
	// Step 1
	// std::cout << "dt_ " << dt_ << std::endl;	
	// std::cout << "Rij " << Rij << std::endl;
	// std::cout << "Aij " << Aij << std::endl;	
	std::cout << "Step 1 \n" << std::endl;
	std::cout << "Dij0 " << Dij0 << std::endl;	
	std::cout << "variable0 = " << variable0 << std::endl;


	variable(ic , jc) = variable0 - alpha1 * dtByAij * (Rij - Dij0);
	variable.correctFlux_f( ic , jc);
	variable.correctFlux_g( ic , jc);

	// if (i == 0 || i == Nci_ - 1)
	// {
	// 	correctWall();
	// }
	
	// if (j == Mci_ - 1 )
	// {
	// 	correctOutlet();
	// }
	// else if (j == 0)
	// {
	// 	correctInlet();
	// }



	variable.print() ;
	variable.printFlux_f() ;
	variable.printFlux_g() ;

	Rij = variable.computeResidualij(i , j);
	
	variable.R().print();
	
	// std::cout << "\n";

	// variable.print() ;
	// variable.printFlux_f() ;
	// variable.printFlux_g() ;
	// std::cout << "R "  <<  std::endl;
	// variable.R().print();

	std::cout << "\n";
	// Step 2
	std::cout << "Step 2 \n" << std::endl;
	// std::cout << "Rij " << Rij << std::endl;	
	std::cout << "Dij0 " << Dij0 << std::endl;	
	// std::cout << "Aij " << Aij << std::endl;	
	variable(ic , jc) = variable0 - alpha2 * dtByAij * (Rij - Dij0);
	variable.correctFlux_f(ic , jc);
	variable.correctFlux_g(ic , jc);
	// if (i == 0 || i == Nci_ - 1)
	// {
	// 	correctWall();
	// }
	
	// if (j == Mci_ - 1 )
	// {
	// 	correctOutlet();
	// }
	// else if (j == 0)
	// {
	// 	correctInlet();
	// }
	variable.print() ;
	variable.printFlux_f() ;
	variable.printFlux_g() ;
	
	Rij = variable.computeResidualij(i , j);
	variable.R().print();
	//std::cout << "\n";

	// variable.print() ;
	// variable.printFlux_f() ;
	// variable.printFlux_g() ;
	// std::cout << "R "  <<  std::endl;
	// variable.R().print();

	std::cout << "\n";

	// Step 3
	std::cout << "Step 3 \n" << std::endl;
	// std::cout << "Rij " << Rij << std::endl;	
	// std::cout << "Dij0 " << Dij0 << std::endl;	
	// std::cout << "Aij " << Aij << std::endl;	
	variable(ic , jc) = variable0 - alpha3 * dtByAij * (Rij - Dij0);
	variable.correctFlux_f(ic , jc);
	variable.correctFlux_g(ic , jc);
	// if (i == 0 || i == Nci_ - 1)
	// {
	// 	correctWall();
	// }
	
	// if (j == Mci_ - 1 )
	// {
	// 	correctOutlet();
	// }
	// else if (j == 0)
	// {
	// 	correctInlet();
	// }
	variable.print() ;
	variable.printFlux_f() ;
	variable.printFlux_g() ;

	Rij = variable.computeResidualij(i , j);
	variable.R().print();
	std::cout << "\n";

	// variable.print() ;
	// variable.printFlux_f() ;
	// variable.printFlux_g() ;
	//std::cout << "R "  <<  std::endl;
	//variable.R().print();

	//std::cout << "\n";

	// Step 4
	std::cout << "Step 4 \n" << std::endl;
	
	variable(ic , jc) = variable0 - alpha4 * dtByAij * (Rij - Dij0);
	
	variable.print() ;
	variable.printFlux_f() ;
	variable.printFlux_g() ;

	
	std::cout << "\n " << std::endl;

}


void Problem::correctTimeStep()
{
	double minDt = 1;

	for (int i = 0; i < Nci_; ++i)
	{
		for (int j = 0; j < Mci_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;

			double lambdaW = 0.5 * (U_(ic , jc) + U_(ic , jc - 1)) + c_(ic , jc);
			double lambdaE = 0.5 * (U_(ic , jc) + U_(ic , jc + 1)) + c_(ic , jc);
			double lambdaN = 0.5 * (V_(ic - 1 , jc) + V_(ic , jc  )) + c_(ic , jc);
			double lambdaS = 0.5 * (V_(ic + 1 , jc) + V_(ic , jc  )) + c_(ic , jc);
			double sumLambda = lambdaW + lambdaE + lambdaN + lambdaS;
			double minDtij = 2 * mesh_.area_(ic , jc)/(sumLambda);

			minDt = std::min(minDt , minDtij);
		}
	}

	dt_ = CFL_ * minDt;

	
}

void Problem::correctInlet()
{
	if (DEBUG)
	{
	}

	std::cout << "Correcting inlet boundary conditions " << std::endl;
	
	double Vinf = Minf_ * cInf_;
	double Riem1_inf = Vinf + 2/(gamma_ - 1) * cInf_;

	for (unsigned int ic = 1; ic < Nc_ - 1; ++ic)
	{
		unsigned int jb = 1;
		// Second option

		// First option
		// p_(ic , 1) = pInf_;
		
		Riem1_(ic , jb) = Riem1_inf;
		Riem2_(ic , jb) = Umag_(ic , jb + 1) - 2/(gamma_- 1)*c_(ic , jb + 1); // From internal domain
		

		Umag_(ic , jb) = 0.5 * (Riem1_(ic , jb) + Riem2_(ic , jb));
		U_(ic , jb) = Umag_(ic , jb);
		c_(ic , jb) = 0.25 * (gamma_ - 1)*(Riem1_(ic, jb) - Riem2_(ic , jb));
		M_(ic , jb) = Umag_(ic , jb)/c_(ic, jb);
		p_(ic , jb) = pInf_/pow((1 + 0.5*(gamma_ - 1) * M_(ic,jb) * M_(ic,jb)), gamma_/(gamma_ - 1));
		rho_(ic , jb) = gamma_ * p_(ic , jb)/c_(ic , jb)/c_(ic , jb);

		// Correct state vector
		rhoU_(ic , jb) = rho_(ic , jb)*U_(ic , jb);
		rhoV_(ic , jb) = rho_(ic , jb)*V_(ic , jb);
		rhoE_(ic , jb) = p_(ic , jb)/(gamma_ - 1) + 0.5 * rho_(ic , jb) * Umag_(ic , jb)*Umag_(ic , jb);
	}

}

void Problem::correctWall()
{
	if (DEBUG)
	{
		std::cout << "Correcting wall boundary conditions " << std::endl;
	}

	for (unsigned int j = 0; j < Mci_; ++j)
	{
		unsigned int  jc = j + 2;
		unsigned int& Icmax = mesh_.Icmax_;
		std::cout << "Icmax " << Icmax << std::endl;

		rho_(Icmax - 1 , jc)  =   rho_(Icmax - 2, jc);	
		rhoE_(Icmax - 1 , jc) =   rhoE_(Icmax - 2, jc);	
		U_(Icmax - 1 , jc)    =   U_(Icmax - 2, jc);	
		V_(Icmax - 1 , jc)    = - V_(Icmax - 2, jc);	
		p_(Icmax - 1 , jc)    =   p_(Icmax - 2, jc);
		c_(Icmax - 1 , jc)    =   c_(Icmax - 2, jc);

		rhoU_(Icmax - 1 , jc) = rho_(Icmax - 1 , jc) * U_(Icmax - 1 , jc);	
		rhoV_(Icmax - 1 , jc) = rho_(Icmax - 1 , jc) * V_(Icmax - 1 , jc);	

		rho_(1 , jc)  =   rho_(2 , jc);	
		rhoE_(1 , jc) =   rhoE_(2 , jc);	
		U_(1 , jc)    =   U_(2 , jc);	
		V_(1 , jc)    = - V_(2 , jc);	
		p_(1 , jc)    =   p_(2 , jc);
		c_(1 , jc)    =   c_(2 , jc);
		rhoU_(1 , jc) = rho_(1 , jc) * U_(1 , jc);	
		rhoV_(1 , jc) = rho_(1 , jc) * V_(1 , jc);	
	}

}

// void Problem::correctPeriodic()
// {
	
// }

void Problem::correctOutlet()
{
	if (DEBUG)
	{
		std::cout << "Correcting outlet boundary conditions " << std::endl;
	}
		// unsigned int ic = i + 2;
		unsigned int Jcmax = Mc_ - 1; 
		std::cout << "Jcmax " << Mc_ - 1 << std::endl;
	for (unsigned int ic = 1; ic < Nc_ - 1; ++ic)
	{
		

		// Zero gradient
		p_(ic , Jcmax - 1) = pInf_*pRatio;
		p_(ic , Jcmax    ) = pInf_*pRatio;

		//Zero gradient
		rho_(ic , Jcmax - 1) = 2 * rho_(ic , Jcmax - 2) - rho_(ic , Jcmax - 3);
		rho_(ic , Jcmax    ) =     rho_(ic , Jcmax - 1);

		rhoU_(ic , Jcmax - 1) = 2 * rhoU_(ic , Jcmax - 2) - rhoU_(ic , Jcmax - 3);
		rhoU_(ic , Jcmax    ) =     rhoU_(ic , Jcmax - 1);

		rhoV_(ic , Jcmax - 1) = 2 * rhoV_(ic , Jcmax - 2) - rhoV_(ic , Jcmax - 3);
		rhoV_(ic , Jcmax    ) =     rhoV_(ic , Jcmax - 1);

		V_(ic , Jcmax - 1) = rhoV_(ic , Jcmax - 1)/rho_(ic , Jcmax - 1);
		V_(ic , Jcmax    ) = V_(ic , Jcmax - 1);

		U_(ic , Jcmax - 1) = rhoU_(ic , Jcmax - 1)/rho_(ic , Jcmax - 1);
		U_(ic , Jcmax    ) = U_(ic , Jcmax - 1);

		Umag_(ic , Jcmax - 1) = sqrt(U_(ic , Jcmax - 1) * U_(ic , Jcmax - 1) + V_(ic , Jcmax - 1) * V_(ic , Jcmax - 1));

		rhoE_(ic , Jcmax - 1) = p_(ic , Jcmax - 1)/(gamma_ - 1) + 0.5 * rho_(ic , Jcmax  - 1) * Umag_(ic , Jcmax - 1) * Umag_(ic , Jcmax - 1);
		rhoE_(ic , Jcmax    ) = rhoE_(ic , Jcmax - 1);


		c_(ic , Jcmax - 1) = sqrt(gamma_ * p_(ic , Jcmax - 1) / rho_(ic , Jcmax - 1));
		c_(ic , Jcmax    ) = sqrt(gamma_ * p_(ic , Jcmax    ) / rho_(ic , Jcmax    ));
		


		// p_(ic , Jcmax - 2 ) = pInf_*pRatio;
		// V_(ic , Jcmax - 2 ) = V_(ic , Jcmax - 3);
		// Riem1_(ic , Jcmax - 3) = sqrt(U_(ic , Jcmax - 3)*U_(ic , Jcmax - 3) + V_(ic , Jcmax - 3)*V_(ic , Jcmax - 3));
		// Riem1_(ic , Jcmax - 2) = Riem1_(ic , Jcmax - 3);

		// double sJcmax_m1 = p_(ic , Jcmax - 3) / pow(rho_(ic , Jcmax - 3) , gamma_);
		// rho_(ic , Jcmax - 2) = pow((p_(ic , Jcmax - 2)/sJcmax_m1) , - gamma_);
		// c_(ic , Jcmax - 2) = sqrt(gamma_ * p_(ic , Jcmax - 2) / rho_(ic , Jcmax - 2) );

		// Umag_(ic , Jcmax - 2) = Riem1_(ic , Jcmax - 2) - 2 * c_(ic , Jcmax - 2)/(gamma_ - 1);

		// U_(ic , Jcmax - 2) =sqrt(Umag_(ic , Jcmax - 2)*Umag_(ic , Jcmax - 2) - V_(ic , Jcmax - 2) * V_(ic , Jcmax - 2));


		// rhoU_(ic , Jcmax - 2) = rho_(ic , Jcmax - 2) * U_(ic , Jcmax - 2); 
		// rhoV_(ic , Jcmax - 2) = rho_(ic , Jcmax - 2) * V_(ic , Jcmax - 2); 
		// rhoE_(ic , Jcmax - 2) = p_(ic , Jcmax - 2)/(gamma_ - 1) + 0.5 * rho_(ic , Jcmax - 2) * Umag_(ic , Jcmax - 2) * Umag_(ic , Jcmax - 2); 







	}	

	// std::cout << "p_(ic , Jcmax - 1)/(gamma_ - 1)" << p_(ic , Jcmax - 1)/(gamma_ - 1) << std::endl;
	// std::cout << "0.5 * rho_(ic , Jcmax  -1) * (U_(ic , Jcmax - 1) * U_(ic , Jcmax - 1) + V_(ic , Jcmax - 1) * V_(ic , Jcmax - 1))" << p_(ic , Jcmax - 1)/(gamma_ - 1) << std::endl;


}

void Problem::correctPropertiesij(unsigned int i , unsigned int j)
{
	U_(i , j) = rhoU_(i , j)/rho_(i , j);
	V_(i , j) = rhoV_(i , j)/rho_(i , j);
	Umag_(i , j) = sqrt(U_(i , j)*U_(i , j) + V_(i , j)*V_(i , j));
	p_(i , j) = (gamma_ - 1) * (rhoE_(i , j) - 0.5 * rho_(i , j) * Umag_(i , j) * Umag_(i , j) );
	c_(i , j) = sqrt( gamma_ * p_(i , j) * rho_(i , j));
}

void Problem::correctProperties()
{
	for (int i = 0; i < Nci_ ; ++i)
	{
		for (int j = 0; j < Mci_ ; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;

			U_(ic , jc) = rhoU_(ic , jc)/rho_(ic , jc);
			V_(ic , jc) = rhoV_(ic , jc)/rho_(ic , jc);
			Umag_(ic , jc) = sqrt(U_(ic , jc)*U_(ic , jc) + V_(ic , jc)*V_(ic , jc));
			p_(ic , jc) = (gamma_ - 1) * (rhoE_(ic , jc) - 0.5 * rho_(ic , jc) * Umag_(ic , jc) * Umag_(ic , jc) );
			c_(ic , jc) = sqrt( gamma_ * p_(ic , jc) * rho_(ic , jc));
			M_(ic , jc) = Umag_(ic , jc)/c_(ic , jc);
		}
	}
}



