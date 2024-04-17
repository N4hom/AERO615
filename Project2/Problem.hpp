#include "Matrix.hpp"
#include "Mesh.hpp"
#include "Variable.hpp"
#include "debug.hpp"
#include <math.h>

class Problem
{

public:
	Mesh mesh_;
	double dt_ = 1e-3;
	unsigned int Nci_ , Mci_;
	unsigned int Nc_ , Mc_;

	Matrix<double> U_, V_ , p_, T_ , Umag_;							   // Primitive variables 
	Variable        rho_, rhoU_,rhoV_, rhoE_;  // State vector components
	Matrix<double> c_;
	Matrix<double> M_;
	Matrix<double> Riem1_;
	Matrix<double> Riem2_;

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
	void computeFluxes();
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
U_(Nc_ , Mc_),
V_(Nc_ , Mc_),
p_(Nc_ , Mc_),
T_(Nc_ , Mc_),
c_(Nc_ , Mc_),
M_(Nc_ , Mc_),
Umag_(Nc_ , Mc_),
Riem1_(Nc_ , Mc_),
Riem2_(Nc_ , Mc_),
rho_("rho"   , Nc_ , Mc_, mesh_, U_, V_, c_ , p_),  // N - 1 being the number of internal cells
rhoU_("rhoU" , Nc_ , Mc_, mesh_, U_, V_, c_ , p_),
rhoV_("rhoV" , Nc_ , Mc_, mesh_, U_, V_, c_ , p_),
rhoE_("rhoE" , Nc_ , Mc_, mesh_, U_, V_, c_ , p_)
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

	computeFluxes();
	

	if(DEBUG)
	{	
		rho_.printFlux_f();

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

	std::cout << "Initializing internal field " << std::endl;
	// Initialize internal field
	for (unsigned int i = 0; i < Nci_; ++i)
	{
		for (unsigned int j = 0; j < Mci_; ++j)
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
			rhoE_(ic , jc) = p_(ic , jc)/(gamma_ - 1) + 0.5 * rho_(ic , jc) * Umag_(ic , jc)*Umag_(ic , jc);
		}
	}

	std::cout << "Initializing inlet " << std::endl;
	//Initialize inlet boundary
	for (unsigned int ic = 1; ic < Nc_ - 1; ++ic)
	{
		for (int j = 0 ; j < 3; ++j)
		{
				
			// unsigned ic = i + 2; 
			// I don't use jc here because there's no need to since I'm setting the boundary conditions to
			// the initial values. Only the first two columns of cells are touched.

			rho_(ic , j)  = rhoInf_;
			p_(ic , j)    = pInf_;
			c_(ic , j)    = cInf_;
			U_(ic , j )   = Minf_ * cInf_;
			rhoU_(ic , j) = rho_(ic , j) * U_(ic , j);
			rhoV_(ic , j) = rho_(ic , j) * V_(ic , j);
			Umag_(ic , j) = sqrt(U_(ic , j)*U_(ic , j) + V_(ic , j) * V_(ic , j));
			rhoE_(ic , j) = p_(ic , j) / (gamma_ - 1) + 0.5 * rho_(ic , j) * Umag_(ic , j)*Umag_(ic , j);
		}
		
	}
	correctInlet();
	correctWall();
	correctOutlet();
}

void Problem::computeFluxes()
{
	if (DEBUG)
	{
		std::cout << "Computing fluxes " << std::endl;
	}
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


	rho_.printFlux_f();

}


void Problem::solve()
{

	unsigned int N = 30;
	unsigned int iter = 0;


	while(iter < N)
	{
		std::cout << "iter " << iter << std::endl;
		rho_.printFlux_f();
		rho_.computeResidual();
		//rho_.computeDissipation();
		rhoU_.computeResidual();
		//rhoU_.computeDissipation();
		rhoV_.computeResidual();
		//rhoV_.computeDissipation();
		rhoE_.computeResidual();
		//rhoE_.computeDissipation();

		double deltaT_ = 1e-5;
		for (int i = 0; i < Nci_; ++i)
		{
			for (int j = 0; j < Mci_; ++j)
			{
				RungeKutta(rho_  , i , j);

				std::cout << "\n " << std::endl;
				//double& Aij       = area(ic,jc);
			}
		}

		std::cout << "After Runge-Kutta for rho " << std::endl;
		rho_.print();

		for (int i = 0; i < Nci_; ++i)
		{
			for (int j = 0; j < Mci_; ++j)
			{
				RungeKutta(rhoU_  , i , j);

				std::cout << "\n " << std::endl;
				//double& Aij       = area(ic,jc);
			}
		}

		for (int i = 0; i < Nci_; ++i)
		{
			for (int j = 0; j < Mci_; ++j)
			{
				RungeKutta(rhoV_  , i , j);

				std::cout << "\n " << std::endl;
				//double& Aij       = area(ic,jc);
			}
		}

		for (int i = 0; i < Nci_; ++i)
		{
			for (int j = 0; j < Mci_; ++j)
			{
				RungeKutta(rhoE_  , i , j);

				std::cout << "\n " << std::endl;
				//double& Aij       = area(ic,jc);
			}
		}
				
		

		
		computeFluxes();

		correctProperties();
		
		correctWall();
		correctInlet();
		correctOutlet();

		std::cout << "End of time step " << std::endl;
		if (DEBUG)
		{
			rho_.print();
			rho_.printFlux_f();
			rho_.printFlux_g();
			rhoU_.print();
			// std::cout << "rhoV " << std::endl;
			rhoV_.print();
			// std::cout << "rhoE " << std::endl;
			rhoE_.print();

			// std::cout << "c " << std::endl;
			// c_.print();

			// std::cout << "M " << std::endl;
			// M_.print();

			std::cout << "Error for rho " << rho_.computeError() << std::endl;
		}

		iter++;
	}





}


void Problem::RungeKutta(Variable& variable, int i , int j )
{
	std::cout << "\n " << std::endl;
	std::cout << "Runge-Kutta " << "variable " << variable.name() << std::endl;
	std::cout << "i " << i << std::endl;
	std::cout << "j " << j << std::endl;

	double alpha1 = 0.25, alpha2 = 0.5 , alpha3 = 1/3, alpha4 = 1;

	Matrix<double>& area = mesh_.area_;

	int ic = i + 2;
	int jc = j + 2;

	double& Aij       = area(ic,jc);
	double& Rij       = variable.R()(i , j);
	double  Dij0      = variable.D()(i , j);
	double  variable0 = variable(ic , jc);

	std::cout << "R "  <<  std::endl;
	variable.R().print();
	// At the next RK step the fluxes must be updated to calculate the new residual. It will be changed only the flux and residual at the cell (ic jc) or (i j)
	
	// Step 1
	std::cout << "Step 1 \n" << std::endl;
	std::cout << "variable0 = " << variable0 << std::endl;
	std::cout << "dt_ " << dt_ << std::endl;	
	std::cout << "Rij " << Rij << std::endl;	
	std::cout << "Dij0 " << Dij0 << std::endl;	
	std::cout << "Aij " << Aij << std::endl;	
	variable(ic , jc) = variable0 - alpha1 * dt_/Aij * (Rij - Dij0);
	std::cout << "updated variable " << variable(ic , jc)  << std::endl;

	variable.correctFlux_f( ic , jc);
	variable.correctFlux_g( ic , jc);
	Rij = variable.computeResidualij(i , j);

	std::cout << "\n";
	// Step 2
	std::cout << "Step 2 \n" << std::endl;
	std::cout << "Rij " << Rij << std::endl;	
	std::cout << "Dij0 " << Dij0 << std::endl;	
	std::cout << "Aij " << Aij << std::endl;	
	variable(ic , jc) = variable0 - alpha2 * dt_/Aij * (Rij - Dij0);
	std::cout << "updated variable " << variable(ic , jc) << std::endl;

	variable.correctFlux_f(ic , jc);
	variable.correctFlux_g(ic , jc);
	Rij = variable.computeResidualij(i , j);


	// Step 3
	std::cout << "Step 3 \n" << std::endl;
	std::cout << "Rij " << Rij << std::endl;	
	std::cout << "Dij0 " << Dij0 << std::endl;	
	std::cout << "Aij " << Aij << std::endl;	
	variable(ic , jc) = variable0 - alpha3 * dt_/Aij * (Rij - Dij0);
	std::cout << "updated variable " << variable(ic , jc)  << std::endl;
	variable.correctFlux_f(ic , jc);
	variable.correctFlux_g(ic , jc);
	Rij = variable.computeResidualij(i , j);


	// Step 4
	std::cout << "Step 4 \n" << std::endl;
	std::cout << "Rij " << Rij << std::endl;	
	std::cout << "Dij0 " << Dij0 << std::endl;	
	std::cout << "Aij " << Aij << std::endl;	
	variable(ic , jc) = variable0 - alpha4 * dt_/Aij * (Rij - Dij0);
	std::cout << "updated variable " << variable(ic , jc)  << std::endl;
	// variable.correctFlux_f(ic , jc);
	// variable.correctFlux_g(ic , jc);
	// Rij = variable.computeResidualij(i , j);

	std::cout << "\n " << std::endl;

}




void Problem::correctInlet()
{
	if (DEBUG)
	{
		std::cout << "Correcting inlet boundary conditions " << std::endl;
	}

	double Vinf = Minf_ * cInf_;
	double Riem1_inf = Vinf + 2/(gamma_-1)*cInf_;

	for (unsigned int ic = 1; ic < Nc_-1; ++ic)
	{
		// unsigned int ic = i + 2;
		
		Riem1_(ic , 2) = Riem1_inf;
		Riem2_(ic , 2) = Umag_(ic , 3) - 2/(gamma_- 1)*c_(ic , 3);
		

		Umag_(ic , 2) = 0.5 * (Riem1_(ic , 2) + Riem2_(ic , 2));
		U_(ic , 2) = Umag_(ic , 2);
		c_(ic , 2) = 0.25 * (gamma_ - 1)*(Riem1_(ic, 2) - Riem2_(ic , 2));
		M_(ic , 2) = Umag_(ic , 2)/c_(ic, 2);
		p_(ic , 2) = pInf_/pow((1 + 0.5*(gamma_ - 1) * M_(ic,2) * M_(ic,2)), gamma_/(gamma_ - 1));
		rho_(ic , 2) = gamma_ * p_(ic , 2)/c_(ic , 1)/c_(ic , 2);
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
		unsigned int Jcmax = mesh_.Jcmax_; // First ghost cell
	for (unsigned int ic = 1; ic < Nc_ - 1; ++ic)
	{

		p_(ic , Jcmax - 1) = pInf_;   // if subsonic!!!!  
		p_(ic , Jcmax    ) = pInf_;
		
		rho_(ic , Jcmax - 1) = 2*rho_(ic , Jcmax - 2) - rho_(ic , Jcmax - 3);
		rho_(ic , Jcmax    ) = rho_(ic , Jcmax - 1);

		c_(ic , Jcmax - 1) = sqrt(gamma_ * p_(ic , Jcmax - 1) / rho_(ic , Jcmax - 1));
		c_(ic , Jcmax    ) = sqrt(gamma_ * p_(ic , Jcmax    ) / rho_(ic , Jcmax    ));

		//V_(ic , Jmax - 1) = V_(ic , Jmax ); 

		rhoU_(ic , Jcmax - 1) = 2*rhoU_(ic , Jcmax - 2) - rhoU_(ic , Jcmax - 3);
		rhoU_(ic , Jcmax    ) = rhoU_(ic , Jcmax - 1);

		rhoV_(ic , Jcmax - 1) = 2*rhoV_(ic , Jcmax - 2) - rhoV_(ic , Jcmax - 3);
		rhoV_(ic , Jcmax    ) = rhoV_(ic , Jcmax - 1);

		V_(ic , Jcmax - 1) = rhoV_(ic , Jcmax - 1)/rho_(ic , Jcmax - 1);
		V_(ic , Jcmax    ) = V_(ic , Jcmax - 1);

		U_(ic , Jcmax - 1) = rhoU_(ic , Jcmax - 1)/rho_(ic , Jcmax - 1);
		U_(ic , Jcmax    ) = U_(ic , Jcmax - 1);


		rhoE_(ic , Jcmax - 1) = p_(ic , Jcmax - 1)/(gamma_ - 1) + 0.5 * rho_(ic , Jcmax  -1) * (U_(ic , Jcmax - 1) * U_(ic , Jcmax - 1) + V_(ic , Jcmax - 1) * V_(ic , Jcmax - 1));
		rhoE_(ic , Jcmax    ) = rhoE_(ic , Jcmax - 1);
	}	

	// std::cout << "p_(ic , Jcmax - 1)/(gamma_ - 1)" << p_(ic , Jcmax - 1)/(gamma_ - 1) << std::endl;
	// std::cout << "0.5 * rho_(ic , Jcmax  -1) * (U_(ic , Jcmax - 1) * U_(ic , Jcmax - 1) + V_(ic , Jcmax - 1) * V_(ic , Jcmax - 1))" << p_(ic , Jcmax - 1)/(gamma_ - 1) << std::endl;


}

void Problem::correctProperties()
{
	for (int i = 0; i < Nci_; ++i)
	{
		for (int j = 0; j < Mci_; ++j)
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

