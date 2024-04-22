#include "Matrix.hpp"
#include "Vector.hpp"
#include "debug.hpp"
#include <algorithm>

// // Storage for coefficients for a single CV
// struct Coefficients
// {
//   double p = 0, n = 0, e = 0, s = 0, w = 0, b = 0;
//   void print(const unsigned int pr = 5) const
//   {
//     cout << setprecision(pr) << scientific << "n = " << n << ", e = " << e << ", s = " << s
//          << ", w = " << w << ", p = " << p << ", b = " << b << endl;
//   }
// };

class Variable
{
	std::string 	 name_;
	unsigned int    Nci_;
	unsigned int    Mci_;
	unsigned int    Nc_;
	unsigned int    Mc_;
	Matrix<double>  phi_;
	Matrix<double>  phiPrev_;
	Matrix<double>  flux_f_;
	Matrix<double>  flux_g_;
	Matrix<Coefficients> Rf_;
	Matrix<Coefficients> Rg_;
	Matrix<Coefficients> DCoeff_;
	Matrix<double>  R_;
	Matrix<double>  D_;
	Mesh& 			mesh_;
	Matrix<double>& p_;
	Matrix<double>& c_;
	Matrix<double>& U_; 
	Matrix<double>& V_;
	Matrix<Coefficients>& s2_;
	Matrix<Coefficients>& s4_;
	Matrix<Coefficients> lambda_;
	double totalResidual_;

public:
	Variable(std::string name, unsigned int N , unsigned int M , Mesh& mesh, Matrix<double>& U , Matrix<double>& V , Matrix<double>& c , Matrix<double>& p , Matrix<Coefficients>& s2_ , Matrix<Coefficients>& s4_);
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
	Matrix<double>& phiPrev(){return phiPrev_;}
	Matrix<double>& flux_f(){return flux_f_;}
	Matrix<double>& flux_g(){return flux_g_;}
	
};

// input is the number of total cell. Easier to set the Variable input using the size of the matrix phi
Variable::Variable(std::string name, unsigned int N, unsigned int M , Mesh& mesh , Matrix<double>& U , Matrix<double>& V , Matrix<double>& c , Matrix<double>& p, Matrix<Coefficients>& s2 , Matrix<Coefficients>& s4):
name_(name),
Nci_(N - 4),
Mci_(M - 4),
Nc_(N),
Mc_(M),
phi_(Nc_ , Mc_),
phiPrev_(Nc_ , Mc_),
flux_f_(Nc_ , Mc_),
flux_g_(Nc_ , Mc_),
Rf_(Nci_ , Mci_),
Rg_(Nci_ , Mci_),
DCoeff_(Nci_ , Mci_),
R_(Nci_ , Mci_),
D_(Nci_ , Mci_),
mesh_(mesh),
p_(p),
c_(c),
U_(U),
V_(V),
s2_(s2),
s4_(s4),
lambda_(Nc_ , Mc_),
totalResidual_(0.)
{}

Variable::~Variable()
{
	
}


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
	return 0.5 * (flux(i - 1  , j   ) + flux(i , j)) ;
}

double Variable::interpolateBottom(Matrix<double>& flux , unsigned int i, unsigned int j) const
{
	return 0.5 * (flux(i + 1  , j   ) + flux(i , j)) ;
}

// The flux should consider the contribution of pressure in rhoUU rhoVV and rhoE
double Variable::flux_f(unsigned int i , unsigned int j)
{
	
	
	if (name_ == "rhoU")
	{
		//std::cout << "phi_ * U + p_" << phi_(i , j) << " " << U_(i , j) << " + " << p_(i , j) << " = " << phi_(i , j) * U_(i , j) + p_(i , j) << endl;
		return phi_(i , j) * U_(i , j) + p_(i , j);
	}
	else if (name_ == "rhoE")
	{
		//std::cout << "phi_ * U + p_ * U_" << phi_(i , j) << " " << U_(i , j) << " + " << p_(i , j) << " * " << U_(i , j) << " = " << phi_(i , j) * U_(i , j) + p_(i , j) * U_(i , j) << endl;
		return phi_(i , j) * U_(i , j) + p_(i , j) * U_(i , j);
	}
	else
	{
		
		
		return phi_(i , j) * U_(i , j);	
	}
}

double Variable::flux_g( unsigned int i , unsigned int j)
{
	

	if (name_ == "rhoV")
	{
		// std::cout << "phi_ * V + p_" << phi_(i , j) << " " << V_(i , j) << " + " << p_(i , j) << " = " << phi_(i , j) * V_(i , j) + p_(i , j) << std::endl;;
		return phi_(i , j) * V_(i , j) + p_(i , j);
	}
	else if (name_ == "rhoE")
	{
		// std::cout << "phi_ * V + p_ * V_" << phi_(i , j) << " " << V_(i , j) << " + " << p_(i , j) << " * " << V_(i , j) << " =" << phi_(i , j) * V_(i , j) + p_(i , j) * U_(i , j) << std::endl;
		return phi_(i , j) * V_(i , j) + p_(i , j) * V_(i , j) ;
	}
	else
	{
		// std::cout << "phi_ * V " << phi_(i , j) << " " << V_(i , j) << " = " << phi_(i , j) * V_(i , j) << std::endl;
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
	// if (DEBUG)
	// {
	// 	std::cout << "Correct flux g " << std::endl; 
	// 	std::cout << "phi_ * V " << phi_(i , j) << " " << V_(i , j) << " = " << phi_(i , j) * V_(i , j) << std::endl;
	// 	/* code */
	// }

	flux_g_(i , j) = flux_g(i , j);

}

void Variable::correctFlux_f(unsigned int i, unsigned int j)
{
	
	// if (DEBUG)
	// {
	// 	std::cout << "Correct flux f " << std::endl; 
	// 	/* code */
	// }
	flux_f_(i , j) = flux_f(i , j);

}


double Variable::computeResidualij(unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;

	// if (DEBUG)
	// {
	// 	std::cout << "xFaces_(i , j).w " << mesh_.xFaces_(i , j).w << std::endl;
	// 	/* code */
	// }
	

	Rf_(i , j).w = interpolateLeft(flux_f_ , ic , jc) * mesh_.xFaces_(i , j).w ;
	Rf_(i , j).e = interpolateRight(flux_f_ , ic , jc) * mesh_.xFaces_(i , j).e;
	Rf_(i , j).n = interpolateTop(flux_f_ , ic , jc) * mesh_.xFaces_(i , j).n ;
	Rf_(i , j).s = interpolateBottom(flux_f_ , ic , jc) * mesh_.xFaces_(i , j).s;

	Rg_(i , j).w = interpolateLeft(flux_g_ , ic , jc) * mesh_.yFaces_(i , j).w ;
	Rg_(i , j).e = interpolateRight(flux_g_ , ic , jc) * mesh_.yFaces_(i , j).e;
	Rg_(i , j).n = interpolateTop(flux_g_ , ic , jc) * mesh_.yFaces_(i , j).n ;
	Rg_(i , j).s = interpolateBottom(flux_g_ , ic , jc) * mesh_.yFaces_(i , j).s;

	


	if (DEBUG)
	{
		std::cout << "f residual " << std::endl;
			Rf_(i , j).print();
		std::cout << "g residual " << std::endl;
			Rg_(i , j).print();
		
		/* code */
	}
	// South and west assumed negative contributions
	double Rij =  (Rf_(i , j).n + Rf_(i , j).s + Rf_(i , j).e + Rf_(i , j).w) - (Rg_(i , j).n + Rg_(i , j).s + Rg_(i , j).e + Rg_(i , j).w);

	Rij = 0.5 * std::abs(mesh_.xFaces_(i , j).e) * (flux_f_(ic , jc + 1) - flux_f_(ic , jc - 1));
	return Rij;
	
}

void Variable::computeResidual()
{	

	if (DEBUG)
	{
		std::cout << "Calculating residual for " << name_ << "\n" << endl;
		this->print();
		this->printFlux_f();
		this->printFlux_g();
		std::cout << " \n " << std::endl;
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

	{
			std::cout << "Residual for " << name_ << "\n" << std::endl;
			R_.print();
	}	
}

void Variable::computeDissipation()
{
	if (DEBUG)
	{
		std::cout << "Calculating dissipation for " << name_ << endl;
	}

	// The face area must be corrected 
	Matrix<Coefficients>& xFaces    = mesh_.xFaces_;
	Matrix<Coefficients>& yFaces    = mesh_.yFaces_;

	// Calculate the second order switches
	if (DEBUG)
	{
		std::cout << "Calculating eigenvalues " << name_  << "\n" << std::endl;
	}

	for (unsigned int i = 0; i < Nci_; ++i)
	{

		for (unsigned int j = 0; j < Mci_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;

			lambda_(ic , jc).n = std::abs(interpolateTop(V_ , ic , jc)) + interpolateTop(c_ , ic , jc);
			lambda_(ic , jc).s = std::abs(interpolateBottom(V_ , ic , jc)) + interpolateBottom(c_ , ic , jc);
			lambda_(ic , jc).e = std::abs(interpolateRight(U_ , ic , jc)) + interpolateRight(c_ , ic , jc);
			lambda_(ic , jc).w = std::abs(interpolateLeft(U_ , ic , jc)) + interpolateLeft(c_ , ic , jc);
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

			
			// The signs should be corrected
			DCoeff_(i , j).w =    s2_(ic , jc).w * std::abs(xFaces(i , j).w) * lambda_(ic , jc).w * (phi_(ic , jc) - phi_(ic , jc - 1)) -
									    s4_(ic , jc).w * std::abs(xFaces(i , j).w) * lambda_(ic , jc).w * ( delta2Eta(phi_ , i , j ) - delta2Eta(phi_ , i , j - 1) ) ;
			DCoeff_(i , j).e =    s2_(ic , jc).e * std::abs(xFaces(i , j).e) * lambda_(ic , jc).e * (phi_(ic , jc + 1) - phi_(ic , jc)) -
								       s4_(ic , jc).e * std::abs(xFaces(i , j).e) * lambda_(ic , jc).e * ( delta2Eta(phi_ , i , j + 1) - delta2Eta(phi_ , i , j) );
			
			DCoeff_(i , j).n =    s2_(ic , jc).n * std::abs(yFaces(i , j).n) * lambda_(ic , jc).n * (phi_(ic  , jc) - phi_(ic - 1, jc)) -
								       s4_(ic , jc).n * std::abs(yFaces(i , j).n) * lambda_(ic , jc).n * ( delta2Csi(phi_ , i , j) - delta2Csi(phi_ , i - 1 , j) );
			DCoeff_(i , j).s =    s2_(ic , jc).s * std::abs(yFaces(i , j).s) * lambda_(ic , jc).s * (phi_(ic + 1 , jc) - phi_(ic  , jc)) -
									    s4_(ic , jc).s * std::abs(yFaces(i , j).s) * lambda_(ic , jc).s * ( delta2Csi(phi_ , i + 1, j) - delta2Csi(phi_ , i , j) );



			//D_(i , j) =  DCoeff_(i , j).s - DCoeff_(i , j).n + DCoeff_(i , j).e - DCoeff_(i , j).w;
			
			D_(i , j) =   DCoeff_(i , j).e - DCoeff_(i , j).w;

			
			
			std::cout << "i " << i << std::endl;
			std::cout << "j " << j << std::endl;
			std::cout << "D Coefficients " << std::endl;
			DCoeff_(i , j).print() ;
			
		}
	}

	std::cout << "Dissipation for " << name_ << "\n" << std::endl;
	D_.print(); 

	
}

double Variable::computeError()
{
	double totalResidual(0.);
	double errMax = 0;
	for (unsigned int i = 0; i < Nci_; ++i)
	{
		for (unsigned int j = 1; j < Mci_; ++j)
		{
			 // std::cout << "D_(i,j) " << D_(i,j) << std::endl;
			 // std::cout << "R_(i,j) " << R_(i,j) << std::endl;
			 // std::cout << "D_(i,j) - R_(i,j) " << D_(i,j) + R_(i,j) << std::endl;
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;

			totalResidual = totalResidual + (D_(i,j) - R_(i,j));
			errMax = max( std::abs(phi_(ic,jc) - phiPrev_(ic,jc)) , errMax)/phi_(ic, jc);





			// std::cout << "totalResidual " << totalResidual << std::endl;

		}
	}

	return errMax;
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
	
	return matrix(ic , jc + 1) - 2*matrix(ic , jc) + matrix(ic , jc - 1);
}