#include "Matrix.hpp"
#include "Vector.hpp"
#include "debug.hpp"
#include <algorithm>

// Storage for coefficients for a single CV
struct Coefficients
{
  double p = 0, n = 0, e = 0, s = 0, w = 0, b = 0;
  void print(const unsigned int pr = 5) const
  {
    cout << setprecision(pr) << scientific << "n = " << n << ", e = " << e << ", s = " << s
         << ", w = " << w << ", p = " << p << ", b = " << b << endl;
  }
};

class Variable
{
	std::string 	 name_;
	unsigned int    Nci_;
	unsigned int    Mci_;
	unsigned int    Nc_;
	unsigned int    Mc_;
	Matrix<double>  phi_;
	Matrix<double>  flux_f_;
	Matrix<double>  flux_g_;
	Matrix<Coefficients> Rf_;
	Matrix<Coefficients> Rg_;
	Matrix<Coefficients> Df_;
	Matrix<Coefficients> Dg_;
	Matrix<double>  R_;
	Matrix<double>  D_;
	Mesh& 			mesh_;
	Matrix<double>& p_;
	Matrix<double>& c_;
	Matrix<double>& U_, V_;
	double totalResidual_;

public:
	Variable(std::string name, unsigned int N , unsigned int M , Mesh& mesh, Matrix<double>& U , Matrix<double>& V , Matrix<double>& c , Matrix<double>& p);
	~Variable();


	double interpolateLeft  (Matrix<double>& flux , unsigned int i, unsigned int j) const;
	double interpolateRight (Matrix<double>& flux , unsigned int i, unsigned int j) const;
	double interpolateTop   (Matrix<double>& flux , unsigned int i, unsigned int j) const;
	double interpolateBottom(Matrix<double>& flux , unsigned int i, unsigned int j) const;

	double flux_f(unsigned int i, unsigned int j);
	double flux_g(unsigned int i, unsigned int j);

	void computeFlux_f();
	void computeFlux_g();
	void correctFlux_f(unsigned int i, unsigned int j);
	void correctFlux_g(unsigned int i, unsigned int j);

	double computeResidualij(unsigned int i , unsigned int j                    );
	void computeResidual();
	void computeDissipation();
	double computeError();
	double Rij(unsigned int i, unsigned int j) const{return R_(i,j);};
	double Dij(unsigned int i, unsigned int j) const{return D_(i,j);};
	Matrix<double>& R(){return R_;};
	Matrix<double>& D(){return D_;};

	double deltaCsi(const Matrix<double>& matrix, unsigned int i, unsigned int j);
	double deltaEta(const Matrix<double>& matrix, unsigned int i, unsigned int j);
	double delta2Csi(const Matrix<double>& matrix, unsigned int i, unsigned int j);
	double delta2Eta(const Matrix<double>& matrix, unsigned int i, unsigned int j);

	std::string name(){return name_;};
	// Const access to the ith,jth element  (read only)
	const double &operator()(unsigned int i, unsigned int j) const;

	// Non const access for getting the ith, jth element
	double &operator()(unsigned int i, unsigned int j);

	void print() const{ std::cout << name_ << std::endl; phi_.print();};
	void printFlux_f() const{ std::cout << "f flux for " << name_ << std::endl; flux_f_.print();};
	void printFlux_g() const{ std::cout << "g flux for " << name_ << std::endl; flux_g_.print();};

	Matrix<double>& phi(){return phi_;}
	Matrix<double>& flux_f(){return flux_f_;}
	Matrix<double>& flux_g(){return flux_g_;}
	
};

// input is the number of total cell. Easier to set the Variable input using the size of the matrix phi
Variable::Variable(std::string name, unsigned int N, unsigned int M , Mesh& mesh , Matrix<double>& U , Matrix<double>& V , Matrix<double>& c , Matrix<double>& p):
name_(name),
Nci_(N - 4),
Mci_(M - 4),
Nc_(N),
Mc_(M),
phi_(Nc_ , Mc_),
flux_f_(Nc_ , Mc_),
flux_g_(Nc_ , Mc_),
Rf_(Nc_ , Mc_),
Rg_(Nc_ , Mc_),
Df_(Nc_ , Mc_),
Dg_(Nc_ , Mc_),
R_(Nci_ , Mci_),
D_(Nci_ , Mci_),
mesh_(mesh),
p_(p),
c_(c),
U_(U),
V_(V),
totalResidual_(0.)
{}

Variable::~Variable()
{}


const double& Variable::operator()(unsigned int i, unsigned int j) const
{
   return phi_(i,j);
}

double& Variable::operator()(unsigned int i, unsigned int j)
{
   return phi_(i,j);
}

double Variable::interpolateLeft(Matrix<double>& flux , unsigned int i, unsigned int j) const
{
	return 0.5 * (flux(i     , j - 1) + flux(i , j)) ;
}

double Variable::interpolateRight(Matrix<double>& flux , unsigned int i, unsigned int j) const
{
	return 0.5 * (flux(i     , j + 1) + flux(i , j)) ;
}

double Variable::interpolateTop(Matrix<double>& flux , unsigned int i, unsigned int j)	const
{
	return 0.5 * (flux(i + 1  , j   ) + flux(i , j)) ;
}

double Variable::interpolateBottom(Matrix<double>& flux , unsigned int i, unsigned int j) const
{
	return 0.5 * (flux(i - 1  , j   ) + flux(i , j)) ;
}

// The flux should consider the contribution of pressure in rhoUU rhoVV and rhoE
double Variable::flux_f(unsigned int i , unsigned int j)
{
	

	if (name_ == "rhoU")
	{
		std::cout << "phi_ * U + p_" << phi_(i , j) << " " << U_(i , j) << " + " << p_(i , j) << " = " << phi_(i , j) * U_(i , j) + p_(i , j) << endl;
		return phi_(i , j) * U_(i , j) + p_(i , j);
	}
	else if (name_ == "rhoE")
	{
		std::cout << "phi_ * U + p_ * U_" << phi_(i , j) << " " << U_(i , j) << " + " << p_(i , j) << " * " << U_(i , j) << " = " << phi_(i , j) * U_(i , j) + p_(i , j) * U_(i , j) << endl;
		return phi_(i , j) * U_(i , j) + p_(i , j) * U_(i , j);
	}
	else
	{
		std::cout << "phi_ * U " << phi_(i , j) << " " << U_(i , j) << " = " << phi_(i , j) * U_(i , j) << std::endl;
		return phi_(i , j) * U_(i , j);	
	}
}

double Variable::flux_g( unsigned int i , unsigned int j)
{
	if (name_ == "rhoV")
	{
		std::cout << "phi_ * V + p_" << phi_(i , j) << " " << V_(i , j) << " + " << p_(i , j) << " = " << phi_(i , j) * V_(i , j) + p_(i , j) << std::endl;;
		return phi_(i , j) * V_(i , j) + p_(i , j);
	}
	else if (name_ == "rhoE")
	{
		std::cout << "phi_ * V + p_ * V_" << phi_(i , j) << " " << V_(i , j) << " + " << p_(i , j) << " * " << V_(i , j) << " =" << phi_(i , j) * V_(i , j) + p_(i , j) * U_(i , j) << std::endl;
		return phi_(i , j) * V_(i , j) + p_(i , j) * V_(i , j) ;
	}
	else
	{
		std::cout << "phi_ * V " << phi_(i , j) << " " << V_(i , j) << " = " << phi_(i , j) * V_(i , j) << std::endl;
		return phi_(i , j) * V_(i , j);	
	}
}



void Variable::computeFlux_f()
{
	for (unsigned int i = 0; i < Nc_; ++i)
	{
		for (unsigned int j = 0; j < Mc_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;
			flux_f_(ic , jc) = flux_f(ic , jc);
		}
	}
}

void Variable::computeFlux_g()
{
	for (unsigned int i = 0; i < Nc_; ++i)
	{
		for (unsigned int j = 0; j < Mc_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;
			flux_g_(ic , jc) = flux_g(ic , jc);
		}
	}
}

void Variable::correctFlux_g(unsigned int i, unsigned int j)
{
	
	std::cout << "flux_g(i , j) " << std::endl;
	flux_g_(i , j) = flux_g(i , j);

}

void Variable::correctFlux_f(unsigned int i, unsigned int j)
{
	
	std::cout << "flux_f(ic , jc) " << std::endl;
	flux_f_(i , j) = flux_f(i , j);

}


double Variable::computeResidualij(unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	
	double RfLeft = interpolateLeft(flux_f_ , ic , jc) * mesh_.xFacesLeft_(i , j) * -1;  // f component calculated assuming that yFacesLeft_ stores face area along csi(y)
	double RfRight = interpolateRight(flux_f_ , ic , jc) * mesh_.xFacesRight_(i , j);
	double RfTop = interpolateTop(flux_f_ , ic , jc) * mesh_.xFacesTop_(i , j) * -1;
	double RfBottom = interpolateBottom(flux_f_ , ic , jc) * mesh_.xFacesBottom_(i , j);

	double RgLeft = interpolateLeft(flux_g_ , ic , jc) * mesh_.yFacesLeft_(i , j) * -1;
	double RgRight = interpolateRight(flux_g_ , ic , jc) * mesh_.yFacesRight_(i , j);
	double RgTop = interpolateTop(flux_g_ , ic , jc) * mesh_.yFacesTop_(i , j) * -1;
	double RgBottom = interpolateBottom(flux_g_ , ic , jc) * mesh_.yFacesBottom_(i , j);


	if (name_ == "rho")
	{
		std::cout << "RfLeft " << RfLeft << std::endl;
		std::cout << "RfRight " << RfRight << std::endl;
		std::cout << "RfTop " << RfTop << std::endl;
		std::cout << "RfBottom " << RfBottom << std::endl;

		std::cout << "RgLeft " << RgLeft << std::endl;
		std::cout << "RgRight " << RgRight << std::endl;
		std::cout << "RgTop " << RgTop << std::endl;
		std::cout << "RgBottom " << RgBottom << std::endl;
	}


	return (RfLeft + RfRight + RfTop + RfBottom) - (RgLeft + RgRight + RgTop + RgBottom);
	
}

void Variable::computeResidual()
{	

	if (DEBUG)
	{
		std::cout << "Calculating residual for " << name_ << "\n" << endl;
	}
	// Residual has the same size as the internal domain, so the loop is over Nci_,Mci_
	for (unsigned int i = 0; i < Nci_; ++i)
	{
		for (unsigned int j = 0; j < Mci_; ++j)
		{
			std::cout << "i " << i << std::endl;
			std::cout << "j " << j << std::endl;
			R_(i,j) = computeResidualij(i,j);
			std::cout << "R_(i,j) = " << R_(i,j) << "\n"<< std::endl;
		}
	}

	
	std::cout << "Residual for " << name_ << "\n" << std::endl;
	R_.print();
	
}

void Variable::computeDissipation()
{
	if (DEBUG)
	{
		std::cout << "Calculating dissipation for " << name_ << endl;
	}
	// const unsigned int Nc = mesh_.Nc_;
	// const unsigned int Mc = mesh_.Mc_;

	const double nu2 = 0;
	const double nu4 = 0.001;

	// The face area must be corrected 
	Matrix<double> yFacesTop    = mesh_.yFacesTop_;
	Matrix<double> yFacesBottom = mesh_.yFacesBottom_;

	Matrix<double> xFacesRight  = mesh_.xFacesRight_;
	Matrix<double> xFacesLeft   = mesh_.xFacesLeft_;

	Matrix<double> sCsi2(Nc_ , Mc_ );
	Matrix<double> sEta2(Nc_ , Mc_ );

	Matrix<double> lambdaCsi(Nc_ , Mc_ );
	Matrix<double> lambdaEta(Nc_ , Mc_ );

	
	// Calculate the second order switches
	if (DEBUG)
	{
		std::cout << "Calculating second order switches for " << name_  << "\n" << std::endl;
	}

	for (unsigned int i = 0; i < Nci_; ++i)
	{

		for (unsigned int j = 0; j < Mci_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;


			// I should be calculating the switch separately at i,j and i+1,j and then take the average
			double sCsi2_ij_num , sCsi2_ij_den;
			double sEta2_ij_num , sEta2_ij_den;

			// delta2csi already adjust the indexes to access the fields
			sCsi2_ij_num = delta2Csi(p_ , i , j);
			sCsi2_ij_den = p_(ic + 1, jc) + 2*p_(ic , jc) + p_(ic - 1, jc);
			double sCsi2_ij = nu2 * sCsi2_ij_num/sCsi2_ij_den;


			sEta2_ij_num = delta2Eta(p_ , i , j);
			sEta2_ij_den = p_(ic , jc + 1) + 2*p_(ic , jc) + p_(ic , jc - 1);
			double sEta2_ij = nu2 * sEta2_ij_num/sEta2_ij_den;

			sCsi2(ic , jc) = sCsi2_ij; 
			sEta2(ic , jc) = sEta2_ij; 


			lambdaCsi(ic , jc) = std::abs(V_(ic , jc)) + c_(ic , jc);
			lambdaEta(ic , jc) = std::abs(U_(ic , jc)) + c_(ic , jc);
		}
	}


	if (DEBUG)
	{
		std::cout << "Calculating fourth order switches and dissipation term "  << name_<< "\n" << std::endl;
	}

	for (unsigned int i = 0; i < Nci_; ++i)
	{
		for (unsigned int j = 0; j < Mci_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;

			double sCsi2Top        = interpolateTop(sCsi2 , ic , jc );
			double sCsi2Bottom     = interpolateBottom(sCsi2 , ic , jc );
			double sCsi2Right      = interpolateRight(sCsi2 , ic , jc );
			double sCsi2Left       = interpolateLeft(sCsi2 , ic , jc );

			double sEta2Top        = interpolateTop(sEta2 , ic , jc );
			double sEta2Bottom     = interpolateBottom(sEta2 , ic , jc );
			double sEta2Right      = interpolateRight(sEta2 , ic , jc );
			double sEta2Left       = interpolateLeft(sEta2 , ic , jc );

			double sCsi4Top        = std::max(0. , nu4 - sCsi2Top);
			double sCsi4Bottom     = std::max(0. , nu4 - sCsi2Bottom);
			double sCsi4Right      = std::max(0. , nu4 - sCsi2Right);
			double sCsi4Left       = std::max(0. , nu4 - sCsi2Left);

			double sEta4Top        = std::max(0. , nu4 - sEta2Top);
			double sEta4Bottom     = std::max(0. , nu4 - sEta2Bottom);
			double sEta4Right      = std::max(0. , nu4 - sEta2Right);
			double sEta4Left       = std::max(0. , nu4 - sEta2Left);

			double lambdaCsiTop    = interpolateTop(lambdaCsi , ic , jc );
			double lambdaCsiBottom = interpolateBottom(lambdaCsi , ic , jc );
			double lambdaCsiRight  = interpolateRight(lambdaCsi , ic , jc );
			double lambdaCsiLeft   = interpolateLeft(lambdaCsi , ic , jc );

			double lambdaEtaTop    = interpolateTop(lambdaEta , ic , jc );
			double lambdaEtaBottom = interpolateBottom(lambdaEta , ic , jc );
			double lambdaEtaRight  = interpolateRight(lambdaEta , ic , jc );
			double lambdaEtaLeft   = interpolateLeft(lambdaEta , ic , jc );

			double deltaCsi2Var = sCsi2Top * yFacesTop(i , j) * lambdaCsiTop * (phi_(ic + 1 , jc) - phi_(ic , jc)) 
								 - sCsi2Bottom * yFacesBottom(i , j) * lambdaCsiBottom * (phi_(ic , jc) - phi_(ic - 1 , jc));
	
			double deltaEta2Var =  sEta2Right * xFacesRight(i , j) * lambdaEtaRight * (phi_(ic , jc + 1) - phi_(ic , jc)) 
							 	 - sEta2Left * xFacesLeft(i , j) * lambdaEtaLeft * (phi_(ic , jc) - phi_(ic , jc - 1));

			double deltaCsi4Var = sCsi4Top * yFacesTop(i , j) * lambdaCsiTop * ( delta2Csi(phi_ , i + 1, j) - delta2Csi(phi_ , i , j) ) 
								 - sCsi4Bottom * yFacesBottom(i , j) * lambdaCsiBottom * ( delta2Eta(phi_ , i , j) - delta2Eta(phi_ , i - 1, j) );

			double deltaEta4Var = sEta4Right * xFacesRight(i , j) * lambdaEtaRight * ( delta2Eta(phi_ , i , j + 1) - delta2Eta(phi_ , i , j) ) 
								 - sEta4Left * xFacesLeft(i , j) * lambdaEtaLeft * ( delta2Eta(phi_ , i , j ) - delta2Eta(phi_ , i , j - 1) );

			D_(i , j)  = (deltaCsi2Var + deltaEta2Var) - (deltaCsi4Var + deltaEta4Var);
			
		
		}
	}

	std::cout << "Dissipation for " << name_ << "\n" << std::endl;
	D_.print(); 

	
}

double Variable::computeError()
{
	double totalResidual(0.);
	for (unsigned int i = 0; i < Nci_; ++i)
	{
		for (unsigned int j = 0; j < Mci_; ++j)
		{
			// std::cout << "D_(i,j) " << D_(i,j) << std::endl;
			// std::cout << "R_(i,j) " << R_(i,j) << std::endl;
			totalResidual = totalResidual + (D_(i,j) - R_(i,j));
			// std::cout << "totalResidual " << totalResidual << std::endl;

		}
	}

	return totalResidual;
}

double Variable::deltaCsi(const Matrix<double>& matrix, unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	double matrixTop = 0.5 * (matrix(ic , jc) + matrix(ic + 1, jc));
	double matrixBottom = 0.5 * (matrix(ic , jc) + matrix(ic - 1, jc));

	return matrixTop - matrixBottom;
}

double Variable::deltaEta(const Matrix<double>& matrix, unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	double matrixRight = 0.5 * (matrix(ic , jc) + matrix(ic , jc + 1));
	double matrixLeft  = 0.5 * (matrix(ic , jc) + matrix(ic , jc - 1));

	return matrixRight - matrixLeft;
}

double Variable::delta2Csi(const Matrix<double>& matrix, unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	
	return matrix(ic + 1, jc) - 2*matrix(ic , jc) + matrix(ic - 1, jc);
}

double Variable::delta2Eta(const Matrix<double>& matrix, unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	
	return matrix(ic , jc + 1) - 2*matrix(ic , jc) + matrix(ic , jc + 1);
}