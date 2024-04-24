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
	double CFL_ = 2;
	const double nu2_ = 0;
	const double nu4_ = 0.1;
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
	double pRatio = 0.99;
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
	void solve();
	void RungeKutta(Variable& variable , unsigned int i , unsigned int j);
	
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
	
	computeFluxes();

	correctInlet();
	correctWall();
	correctOutlet();



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
	double fluxLeft = 0.5 * (flux(i     , j - 1) + flux(i , j));
	return fluxLeft ;
}

double Problem::interpolateRight(Matrix<double>& flux , unsigned int i, unsigned int j) const
{
	double fluxRight = 0.5 * (flux(i     , j + 1) + flux(i , j)) ;
	return fluxRight;
}

double Problem::interpolateTop(Matrix<double>& flux , unsigned int i, unsigned int j)	const
{
	double fluxTop = 0.5 * (flux(i - 1  , j   ) + flux(i , j)) ;
	return fluxTop;
}

double Problem::interpolateBottom(Matrix<double>& flux , unsigned int i, unsigned int j) const
{
	double fluxBottom =  0.5 * (flux(i + 1  , j   ) + flux(i , j)) ;
	return fluxBottom;
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

	
	// Initialize internal field
	for (unsigned int ic = 0; ic < Nc_; ++ic)
	{
		for (unsigned int jc = 0; jc < Mc_ ; ++jc)
		{
			// unsigned ic = i + 2;
			// unsigned jc = j + 2;
			rho_(ic , jc) = rhoInf_;
			p_(ic , jc) = pInf_*pRatio;
			c_(ic , jc) = cInf_;
			U_(ic , jc) = cInf_ * Minf_;
			rhoU_(ic , jc) = rho_(ic , jc) * U_(ic , jc);
			rhoV_(ic , jc) = rho_(ic , jc) * V_(ic , jc);
			Umag_(ic , jc) = sqrt(U_(ic , jc)*U_(ic , jc) + V_(ic , jc) * V_(ic , jc));
			rhoE_(ic , jc) = p_(ic , jc)/(gamma_ - 1) + 0.5 * rho_(ic , jc) * Umag_(ic , jc)*Umag_(ic , jc);
		}
	}

	
}

void Problem::computeSwitches()
{
	for (unsigned int i = - 1; i < Nci_ + 1; ++i)
	{
		for (unsigned int j = - 1; j < Mci_ + 1; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;

			s2Csi_(ic , jc) = nu2_ * std::abs(delta2Csi(p_ , i , j))/(p_(ic , jc + 1) + 2*p_(ic , jc) + p_(ic , jc - 1));
			s2Eta_(ic , jc) = nu2_ * std::abs(delta2Eta(p_ , i , j))/(p_(ic , jc + 1) + 2*p_(ic , jc) + p_(ic , jc - 1));
		
		}
			
	}


	// Interpolate switches at cell faces
	
	for (unsigned int i = 0; i < Nci_ ; ++i)
	{
		for (unsigned int j = 0; j < Mci_; ++j)
		{

			unsigned int ic = i + 2;
			unsigned int jc = j + 2;


			s2_(ic , jc).e = interpolateRight(s2Eta_ , ic , jc);
			s2_(ic  , jc).w = interpolateLeft(s2Eta_ , ic , jc);

			// if ghost cells are used (??)
			s2_(ic , jc).n = interpolateTop(s2Csi_ , ic , jc);
			s2_(ic , jc).s = interpolateBottom(s2Csi_ , ic , jc);

			// if (i == 0)
			// {
			// 	s2_(ic , jc).n = 0;
			// 	s2_(ic , jc).s = interpolateBottom(s2Csi_ , ic , jc);
			// }
			// else if (i == Nci_ - 1)
			// {
			// 	s2_(ic , jc).s = 0;
			// 	s2_(ic , jc).n = interpolateTop(s2Csi_ , ic , jc);

			// }
			// else
			// {
			// 	s2_(ic , jc).n = interpolateTop(s2Csi_ , ic , jc);
			// 	s2_(ic , jc).s = interpolateBottom(s2Csi_ , ic , jc);

			// }

			
			s4_(ic , jc).e = std::max(0. , nu4_ - s2_(ic , jc).e );
			s4_(ic , jc).w = std::max(0. , nu4_ - s2_(ic , jc).w );

			s4_(ic , jc).n = std::max(0. , nu4_ - s2_(ic , jc).n );
			s4_(ic , jc).s = std::max(0. , nu4_ - s2_(ic , jc).s );
			
			// if (i == 0)
			// {
			// 	s4_(ic , jc).n = 0;
			// 	s4_(ic , jc).s = std::max(0. , nu4_ - s2_(ic , jc).s );

			// }
			// else if (i == Nci_ - 1)
			// {

			// 	s4_(ic , jc).n = std::max(0. , nu4_ - s2_(ic , jc).n );
			// 	s4_(ic , jc).s = 0;

			// }
			// else
			// {
			// 	s4_(ic , jc).n = std::max(0. , nu4_ - s2_(ic , jc).n );
			// 	s4_(ic , jc).s = std::max(0. , nu4_ - s2_(ic , jc).s );
			// }


			
		}
			
	}


	

}

void Problem::computeFluxes()
{
	if (DEBUG)
	{
		std::cout << "Computing fluxes " << std::endl;
	}
	
	// Looping over all domain because the fluxes in the ghost cells are needed
	for (unsigned int i = 0; i < Nci_; ++i)
	{
		for (unsigned int j = 0; j < Mci_; ++j)
		{
			unsigned ic = i + 2;
			unsigned jc = j + 2;

			// rho_.correctFlux_f(ic , jc);
			// rho_.correctFlux_g(ic , jc);
			
			// rhoU_.correctFlux_f(ic , jc);
			// rhoU_.correctFlux_g(ic , jc);
			
			// rhoV_.correctFlux_f(ic , jc);
			// rhoV_.correctFlux_g(ic , jc);
			
			// rhoE_.correctFlux_f(ic , jc);
			// rhoE_.correctFlux_g(ic , jc);

			rho_.flux_f()(ic , jc) = rho_(ic , jc) * U_(ic , jc);
			rho_.flux_g()(ic , jc) = rho_(ic , jc) * V_(ic , jc);

			rhoU_.flux_f()(ic , jc) = rhoU_(ic , jc) * U_(ic , jc) + p_(ic , jc);
			rhoU_.flux_g()(ic , jc) = rhoU_(ic , jc) * V_(ic , jc);

			rhoV_.flux_f()(ic , jc) = rhoV_(ic , jc) * U_(ic , jc);
			rhoV_.flux_g()(ic , jc) = rhoV_(ic , jc) * V_(ic , jc) + p_(ic , jc);

			rhoE_.flux_f()(ic , jc) = (rhoE_(ic , jc) + p_(ic , jc)) * U_(ic , jc);
			rhoE_.flux_g()(ic , jc) = (rhoE_(ic , jc) + p_(ic , jc)) * V_(ic , jc) ;

			
		}
	}


	

}


void Problem::solve()
{

	unsigned int N = 100;
	unsigned int iter = 0;


	while(iter < N)
	{
		if (DEBUG)
		{
			std::cout << "-------------------------------------------------------------------" << std::endl;
			std::cout << "iter " << iter << std::endl;
			std::cout << "dt " << dt_ << std::endl;
		}
		
		rho_.print();
		rho_.printFlux_f();
		rho_.printFlux_g();
		rhoU_.print();
		rhoU_.printFlux_f();
		rhoU_.printFlux_g();
		rhoV_.print();
		rhoV_.printFlux_f();
		rhoV_.printFlux_g();
		rhoE_.print();
		rhoE_.printFlux_f();
		rhoE_.printFlux_g();
	

		correctTimeStep();

		computeSwitches();

		rho_.computeResidual();
		rho_.computeDissipation();
		
		rhoU_.computeResidual();
		rhoU_.computeDissipation();
		
		rhoV_.computeResidual();
		rhoV_.computeDissipation();
		
		rhoE_.computeResidual();
		rhoE_.computeDissipation();
		
		for (unsigned int i = 0; i < Nci_; ++i)
		{
			for (unsigned int j = 1; j < Mci_  ; ++j)
			{
				unsigned int ic = i + 2;
				unsigned int jc = j + 2;

				rho_.phiPrev()(ic,jc) = rho_(ic , jc);
				RungeKutta(rho_  , i , j);
				RungeKutta(rhoU_  , i , j);
				RungeKutta(rhoV_  , i , j);
				RungeKutta(rhoE_  , i , j);
				correctProperties();
			}
		}

		if (DEBUG)
		{
			outputFile <<  rho_.computeError() << std::endl;
			
			rho_.print();
			rho_.printFlux_f();
			rho_.printFlux_g();
			rhoU_.print();
			rhoU_.printFlux_f();
			rhoU_.printFlux_g();
			rhoV_.print();
			rhoV_.printFlux_f();
			rhoV_.printFlux_g();
			rhoE_.print();
			rhoE_.printFlux_f();
			rhoE_.printFlux_g();
			std::cout << "p " ;
			p_.print();
			std::cout << "Umag " ;
			Umag_.print();
			std::cout << "u " ;
			U_.print();
			std::cout << "u " ;
			V_.print();
			std::cout << "c " ;
			c_.print();
			std::cout << "Riem1_ " ;
			Riem1_.print();
			std::cout << "Riem2_ " ;
			Riem2_.print();
		}
		
		correctProperties();
		computeFluxes();
		
		std::cout << "Correcting boundary conditions " << std::endl;
		
		
		correctInlet();
		correctWall();
		correctOutlet();

		rho_.print();
		rho_.printFlux_f();
		rho_.printFlux_g();
		rhoU_.print();
		rhoU_.printFlux_f();
		rhoU_.printFlux_g();
		rhoV_.print();
		rhoV_.printFlux_f();
		rhoV_.printFlux_g();
		rhoE_.print();
		rhoE_.printFlux_f();
		rhoE_.printFlux_g();
		
		

		

		std::cout << "End of time step " << std::endl;

		std::cout << "--------------------------------------------------------------------------------" << std::endl;

		iter++;
	}





}


void Problem::RungeKutta(Variable& variable, unsigned int i ,  unsigned int j )
{	
	const double alpha1 = 0.25, alpha2 = 0.5 , alpha3 = 0.33333, alpha4 = 1.;
	
	Matrix<double>& area = mesh_.area_;

	unsigned int ic = i + 2;
	unsigned int jc = j + 2;


	std::cout << "RungeKutta " << variable.name() << std::endl;
	variable.print();
	variable.printFlux_f();
	variable.printFlux_g();
	std::cout << "i " << i << std::endl;
	std::cout << "j " << j << std::endl;
	double& Aij       = area(ic,jc);
	double& Rij       = variable.R()(i , j);
	double  Dij0      = variable.D()(i , j);
	//double  Dij0      = 0;
	double  variable0 = variable(ic , jc);
	double dtByAij = dt_/Aij;
	// At the next RK step the fluxes must be updated to calculate the new residual. It will be changed only the flux and residual at the cell (ic jc) or (i j)
	
	std::cout << "D " << std::endl;
	variable.D().print();	
	std::cout << "R " << std::endl;
	variable.R().print();

	variable(ic , jc) = variable0 - alpha1 * dtByAij * (Rij - Dij0);
	variable.correctFlux_f( ic , jc);
	variable.correctFlux_g( ic , jc);

	Rij = variable.computeResidualij(i , j);
	std::cout << "R " << std::endl;
	variable.R().print();
	
	variable(ic , jc) = variable0 - alpha2 * dtByAij * (Rij - Dij0);
	variable.correctFlux_f(ic , jc);
	variable.correctFlux_g(ic , jc);

	Rij = variable.computeResidualij(i , j);

	std::cout << "R " << std::endl;
	variable.R().print();
	
	variable(ic , jc) = variable0 - alpha3 * dtByAij * (Rij - Dij0);
	variable.correctFlux_f(ic , jc);
	variable.correctFlux_g(ic , jc);
	
	Rij = variable.computeResidualij(i , j);
	
	std::cout << "R " << std::endl;
	variable.R().print();
	
	variable(ic , jc) = variable0 - alpha4 * dtByAij * (Rij - Dij0);
	
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

			double lambdaW = std::abs(0.5 * (U_(ic , jc) + U_(ic , jc - 1))) + c_(ic , jc);
			double lambdaE = std::abs(0.5 * (U_(ic , jc) + U_(ic , jc + 1))) + c_(ic , jc);
			double lambdaN = std::abs(0.5 * (V_(ic - 1 , jc) + V_(ic , jc  ))) + c_(ic , jc);
			double lambdaS = std::abs(0.5 * (V_(ic + 1 , jc) + V_(ic , jc  ))) + c_(ic , jc);
			double sumLambda = lambdaW*mesh_.xFaces_(i , j).w + lambdaE*mesh_.xFaces_(i , j).e + lambdaN*mesh_.yFaces_(i , j).n + lambdaS*mesh_.xFaces_(i , j).s;
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

	for (unsigned int ic = 0 ; ic < Nc_ ; ++ic)
	{
		unsigned int jb = 2;
		// Second option

		// First option
		// p_(ic , 1) = pInf_;
		
		Riem1_(ic , jb) = Riem1_inf;
		Riem2_(ic , jb) = Umag_(ic , jb + 1) - 2/(gamma_- 1)*c_(ic , jb + 1); // From internal domain
		
		// Two layers of ghost cells are needed because of the dissipation scheme
		Umag_(ic , jb) = 0.5 * (Riem1_(ic , jb) + Riem2_(ic , jb));
		Umag_(ic , jb - 1) = Vinf;
		Umag_(ic , jb - 2) = Vinf;
		
		U_(ic , jb) = Umag_(ic , jb);
		U_(ic , jb - 1) = Vinf;
		U_(ic , jb - 2) = Vinf;
		
		c_(ic , jb) = 0.25 * (gamma_ - 1)*(Riem1_(ic, jb) - Riem2_(ic , jb));
		c_(ic , jb - 1) = cInf_;
		c_(ic , jb - 2) = cInf_;
		
		M_(ic , jb) = Umag_(ic , jb)/c_(ic, jb);
		M_(ic , jb - 1) = Minf_;
		M_(ic , jb - 2) = Minf_;
		
		p_(ic , jb) = pInf_/pow((1 + 0.5*(gamma_ - 1) * M_(ic,jb) * M_(ic,jb)), gamma_/(gamma_ - 1));
		p_(ic , jb - 1) = pInf_;
		p_(ic , jb - 2) = pInf_;
		
		rho_(ic , jb) = gamma_ * p_(ic , jb)/(c_(ic , jb)*c_(ic , jb));
		rho_(ic , jb - 1) = rhoInf_;
		rho_(ic , jb - 2) = rhoInf_;

		rho_.correctFlux_f(ic , jb    );
		rho_.correctFlux_f(ic , jb - 1);
		rho_.correctFlux_f(ic , jb - 2);
		rho_.correctFlux_g(ic , jb    );
		rho_.correctFlux_g(ic , jb - 1);
		rho_.correctFlux_g(ic , jb - 2);

		// Correct state vector
		rhoU_(ic , jb) = rho_(ic , jb)*U_(ic , jb);
		rhoU_(ic , jb - 1) = rhoInf_ * Uinf_;
		rhoU_(ic , jb - 2) = rhoInf_ * Uinf_;

		rhoU_.correctFlux_f(ic , jb    );
		rhoU_.correctFlux_f(ic , jb - 1);
		rhoU_.correctFlux_f(ic , jb - 2);
		rhoU_.correctFlux_g(ic , jb    );
		rhoU_.correctFlux_g(ic , jb - 1);
		rhoU_.correctFlux_g(ic , jb - 2);
		
		rhoV_(ic , jb) = rho_(ic , jb)*V_(ic , jb);
		rhoV_(ic , jb - 1) = 0;
		rhoV_(ic , jb - 2) = 0;
	
		rhoV_.correctFlux_f(ic , jb    );
		rhoV_.correctFlux_f(ic , jb  - 1    );
		rhoV_.correctFlux_f(ic , jb  - 2    );
		rhoV_.correctFlux_g(ic , jb    );
		rhoV_.correctFlux_g(ic , jb - 1);
		rhoV_.correctFlux_g(ic , jb - 2);
		
		rhoE_(ic , jb) = p_(ic , jb)/(gamma_ - 1) + 0.5 * rho_(ic , jb) * Umag_(ic , jb)*Umag_(ic , jb);
		rhoE_(ic , jb - 1) = p_(ic , jb - 1)/(gamma_ - 1) + 0.5 * rho_(ic , jb - 1) * Umag_(ic , jb - 1)*Umag_(ic , jb - 1);
		rhoE_(ic , jb - 2) = p_(ic , jb - 2)/(gamma_ - 1) + 0.5 * rho_(ic , jb - 2) * Umag_(ic , jb - 2)*Umag_(ic , jb - 2);

		rhoE_.correctFlux_f(ic , jb     );
		rhoE_.correctFlux_f(ic , jb - 1 );
		rhoE_.correctFlux_f(ic , jb - 2 );
		rhoE_.correctFlux_g(ic , jb     );
		rhoE_.correctFlux_g(ic , jb - 1 );
		rhoE_.correctFlux_g(ic , jb - 2 );

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

		rho_(Icmax - 1 , jc)  =   rho_(Icmax - 2, jc);	
		rho_(Icmax     , jc)  =   rho_(Icmax - 3, jc);

		rhoE_(Icmax - 1 , jc) =   rhoE_(Icmax - 2, jc);	
		rhoE_(Icmax     , jc) =   rhoE_(Icmax - 3, jc);	
		
		U_(Icmax - 1 , jc)    =   U_(Icmax - 2, jc);	
		U_(Icmax     , jc)    =   U_(Icmax - 3, jc);	

		V_(Icmax - 1 , jc)    = - V_(Icmax - 2, jc);	
		V_(Icmax     , jc)    = - V_(Icmax - 3, jc);	

		Umag_(Icmax - 1 , jc)    = sqrt(U_(Icmax - 1 , jc) * U_(Icmax - 1 , jc) + V_(Icmax - 1 , jc) * V_(Icmax - 1 , jc));	
		Umag_(Icmax     , jc)    = sqrt(U_(Icmax     , jc) * U_(Icmax     , jc) + V_(Icmax     , jc) * V_(Icmax     , jc));
	
		p_(Icmax - 1 , jc)    =   p_(Icmax - 2, jc);
		p_(Icmax     , jc)    =   p_(Icmax - 3, jc);


		c_(Icmax - 1 , jc)    =   c_(Icmax - 2, jc);
		c_(Icmax     , jc)    =   c_(Icmax - 3, jc);

		rhoU_(Icmax - 1 , jc) = rho_(Icmax - 1 , jc) * U_(Icmax - 1 , jc);	
		rhoU_(Icmax     , jc) = rho_(Icmax     , jc) * U_(Icmax     , jc);	
	
		rhoV_(Icmax - 1 , jc) = rho_(Icmax - 1 , jc) * V_(Icmax - 1 , jc);	
		rhoV_(Icmax     , jc) = rho_(Icmax     , jc) * V_(Icmax     , jc);	

		rho_(1 , jc)  =   rho_(2 , jc);	
		rho_(0 , jc)  =   rho_(3 , jc);	

		rhoE_(1 , jc) =   rhoE_(2 , jc);	
		rhoE_(0 , jc) =   rhoE_(3 , jc);	


		U_(1 , jc)    =   U_(2 , jc);	
		U_(0 , jc)    =   U_(3 , jc);	

		V_(1 , jc)    = - V_(2 , jc);	
		V_(0 , jc)    = - V_(3 , jc);	

		Umag_(1 , jc)    = sqrt(U_(1 , jc) * U_(1 , jc) + V_(1 , jc) * V_(1 , jc));	
		Umag_(0 , jc)    = sqrt(U_(0 , jc) * U_(0 , jc) + V_(0 , jc) * V_(0 , jc));

		p_(1 , jc)    =   p_(2 , jc);
		p_(0 , jc)    =   p_(3 , jc);
	
		c_(1 , jc)    =   c_(2 , jc);
		c_(0 , jc)    =   c_(3 , jc);
		
		rhoU_(1 , jc) = rho_(1 , jc) * U_(1 , jc);	
		rhoU_(0 , jc) = rho_(0 , jc) * U_(0 , jc);	

		rhoV_(1 , jc) = rho_(1 , jc) * V_(1 , jc);	
		rhoV_(0 , jc) = rho_(0 , jc) * V_(0 , jc);	

		rho_.correctFlux_f(1         , jc);
		rho_.correctFlux_f(0         , jc);
		rho_.correctFlux_g(1         , jc);
		rho_.correctFlux_g(0         , jc);
		rho_.correctFlux_f(Icmax     , jc );
		rho_.correctFlux_f(Icmax - 1 , jc );
		rho_.correctFlux_g(Icmax     , jc );
		rho_.correctFlux_g(Icmax - 1 , jc );

		
		rhoU_.correctFlux_f(1         , jc);
		rhoU_.correctFlux_f(0         , jc);
		rhoU_.correctFlux_g(1         , jc);
		rhoU_.correctFlux_g(0         , jc);
		rhoU_.correctFlux_f(Icmax     , jc );
		rhoU_.correctFlux_f(Icmax - 1 , jc );
		rhoU_.correctFlux_g(Icmax     , jc );
		rhoU_.correctFlux_g(Icmax - 1 , jc );

		rhoV_.correctFlux_f(1         , jc);
		rhoV_.correctFlux_f(0         , jc);
		rhoV_.correctFlux_g(1         , jc);
		rhoV_.correctFlux_g(0         , jc);
		rhoV_.correctFlux_f(Icmax     , jc );
		rhoV_.correctFlux_f(Icmax - 1 , jc );
		rhoV_.correctFlux_g(Icmax     , jc );
		rhoV_.correctFlux_g(Icmax - 1 , jc );

		rhoE_.correctFlux_f(1         , jc);
		rhoE_.correctFlux_f(0         , jc);
		rhoE_.correctFlux_g(1         , jc);
		rhoE_.correctFlux_g(0         , jc);
		rhoE_.correctFlux_f(Icmax     , jc );
		rhoE_.correctFlux_f(Icmax - 1 , jc );
		rhoE_.correctFlux_g(Icmax     , jc );
		rhoE_.correctFlux_g(Icmax - 1 , jc );

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
	
	for (unsigned int ic = 0 ; ic < Nc_ ; ++ic)
	{
		

		// Zero second derivative
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
		Umag_(ic , Jcmax    ) = Umag_(ic , Jcmax - 1);

		rhoE_(ic , Jcmax - 1) = p_(ic , Jcmax - 1)/(gamma_ - 1) + 0.5 * rho_(ic , Jcmax  - 1) * Umag_(ic , Jcmax - 1) * Umag_(ic , Jcmax - 1);
		rhoE_(ic , Jcmax    ) = rhoE_(ic , Jcmax - 1);


		c_(ic , Jcmax - 1) = sqrt(gamma_ * p_(ic , Jcmax - 1) / rho_(ic , Jcmax - 1));
		c_(ic , Jcmax    ) = sqrt(gamma_ * p_(ic , Jcmax    ) / rho_(ic , Jcmax    ));


		rho_.correctFlux_f(ic , Jcmax - 1);
		rho_.correctFlux_f(ic , Jcmax    );
		rho_.correctFlux_g(ic , Jcmax - 1);
		rho_.correctFlux_g(ic , Jcmax    );

		rhoU_.correctFlux_f(ic , Jcmax - 1);
		rhoU_.correctFlux_f(ic , Jcmax    );
		rhoU_.correctFlux_g(ic , Jcmax - 1);
		rhoU_.correctFlux_g(ic , Jcmax    );

		rhoV_.correctFlux_f(ic , Jcmax - 1);
		rhoV_.correctFlux_f(ic , Jcmax    );
		rhoV_.correctFlux_g(ic , Jcmax - 1);
		rhoV_.correctFlux_g(ic , Jcmax    );

		rhoE_.correctFlux_f(ic , Jcmax - 1);
		rhoE_.correctFlux_f(ic , Jcmax    );
		rhoE_.correctFlux_g(ic , Jcmax - 1);
		rhoE_.correctFlux_g(ic , Jcmax    );

	

	}	


}

void Problem::correctPropertiesij(unsigned int i , unsigned int j)
{
	U_(i , j) = rhoU_(i , j)/rho_(i , j);
	V_(i , j) = rhoV_(i , j)/rho_(i , j);
	Umag_(i , j) = sqrt(U_(i , j)*U_(i , j) + V_(i , j)*V_(i , j));
	p_(i , j) = (gamma_ - 1) * (rhoE_(i , j) - 0.5 * rho_(i , j) * Umag_(i , j) * Umag_(i , j) );
	c_(i , j) = sqrt( gamma_ * p_(i , j) / rho_(i , j));
	M_(i , j) = Umag_(i , j)/c_(i , j);
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
			Umag_(ic , jc) = sqrt( U_(ic , jc) * U_(ic , jc) + V_(ic , jc) * V_(ic , jc));
			p_(ic , jc) = (gamma_ - 1) * (rhoE_(ic , jc) - 0.5 * rho_(ic , jc) * Umag_(ic , jc) * Umag_(ic , jc) );
			c_(ic , jc) = sqrt( gamma_ * p_(ic , jc) / rho_(ic , jc));
			M_(ic , jc) = Umag_(ic , jc)/c_(ic , jc);
		}
	}
}



