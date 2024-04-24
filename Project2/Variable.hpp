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
		return phi_(i , j) * U_(i , j) + p_(i , j);
	}
	else if (name_ == "rhoE")
	{
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
		return phi_(i , j) * V_(i , j) + p_(i , j);
	}
	else if (name_ == "rhoE")
	{
		return phi_(i , j) * V_(i , j) + p_(i , j) * V_(i , j) ;
	}
	else
	{
		return phi_(i , j) * V_(i , j);	
	}
}



// void Variable::computeFlux_f()
// {
// 	for (unsigned int i = 0; i < Nci_; ++i)
// 	{
// 		for (unsigned int j = 0; j < Mci_; ++j)
// 		{
// 			unsigned int ic = i + 2;
// 			unsigned int jc = j + 2;
// 			flux_f_(ic , jc) = flux_f(ic , jc);
// 		}
// 	}
// }

// void Variable::computeFlux_g()
// {
// 	for (unsigned int i = 0; i < Nci_; ++i)
// 	{
// 		for (unsigned int j = 0; j < Mci_; ++j)
// 		{
// 			unsigned int ic = i + 2;
// 			unsigned int jc = j + 2;
// 			flux_g_(ic , jc) = flux_g(ic , jc);
// 		}
// 	}
// }

void Variable::correctFlux_g(unsigned int i, unsigned int j)
{
	flux_g_(i , j) = flux_g(i , j);

}

void Variable::correctFlux_f(unsigned int i, unsigned int j)
{
	flux_f_(i , j) = flux_f(i , j);
}


double Variable::computeResidualij(unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;


	
	Rf_(i , j).w = interpolateLeft(flux_f_ , ic , jc) * mesh_.xFaces_(i , j).w;
	Rf_(i , j).e = interpolateRight(flux_f_ , ic , jc) * mesh_.xFaces_(i , j).e;
	Rf_(i , j).n = interpolateTop(flux_f_ , ic , jc) * mesh_.xFaces_(i , j).n ;
	Rf_(i , j).s = interpolateBottom(flux_f_ , ic , jc) * mesh_.xFaces_(i , j).s;

	
	Rg_(i , j).w = interpolateLeft(flux_g_ , ic , jc) * mesh_.yFaces_(i , j).w ;
	Rg_(i , j).e = interpolateRight(flux_g_ , ic , jc) * mesh_.yFaces_(i , j).e;
	Rg_(i , j).n = interpolateTop(flux_g_ , ic , jc) * mesh_.yFaces_(i , j).n ;
	Rg_(i , j).s = interpolateBottom(flux_g_ , ic , jc) * mesh_.yFaces_(i , j).s;

	// South and west assumed negative contributions
	double Rij =  (Rf_(i , j).s - Rf_(i , j).n + Rf_(i , j).e - Rf_(i , j).w) - (Rg_(i , j).s - Rg_(i , j).n + Rg_(i , j).e - Rg_(i , j).w);
	// double Rij =  (Rf_(i , j).s - Rf_(i , j).n + Rf_(i , j).e - Rf_(i , j).w) ;
	if (i == 3 && j == 3)
	{
		std::cout << "Rij " << Rij << std::endl;  
		/* code */
	}

	// Rij = 0.5 * std::abs(mesh_.xFaces_(i , j).e) * (flux_f_(ic , jc + 1) - flux_f_(ic , jc - 1));
	return Rij;
	
}

void Variable::computeResidual()
{	

	// Residual has the same size as the internal domain, so the loop is over Nci_,Mci_
	for (unsigned int i = 0; i < Nci_; ++i)
	{
		for (unsigned int j = 0; j < Mci_; ++j)
		{
			R_(i,j) = computeResidualij(i,j);
		}
	}

}

void Variable::computeDissipation()
{

	// The face area must be corrected 
	Matrix<Coefficients>& xFaces    = mesh_.xFaces_;
	Matrix<Coefficients>& yFaces    = mesh_.yFaces_;


	for (unsigned int i = 0; i < Nci_; ++i)
	{

		for (unsigned int j = 0; j < Mci_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;

			lambda_(ic , jc).n = std::abs(interpolateTop(V_ , ic , jc)) + c_(ic , jc);
			lambda_(ic , jc).s = std::abs(interpolateBottom(V_ , ic , jc)) + c_(ic , jc);
			lambda_(ic , jc).e = std::abs(interpolateRight(U_ , ic , jc)) + c_(ic , jc);
			lambda_(ic , jc).w = std::abs(interpolateLeft(U_ , ic , jc)) + c_ (ic , jc);
		}
	}


	for (unsigned int i = 0; i < Nci_; ++i)
	{
		for (unsigned int j = 0; j < Mci_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;

			double delta2W = phi_(ic , jc) - phi_(ic , jc - 1);
			double delta2E = phi_(ic , jc) - phi_(ic , jc + 1);
			double delta2N = phi_(ic , jc) - phi_(ic - 1, jc);
			double delta2S = phi_(ic + 1, jc) - phi_(ic   , jc);

			double delta4W = phi_(ic     , jc + 1) - 3 * phi_(ic     , jc    ) + 3 * phi_(ic      , jc - 1) - phi_(ic     , jc - 2);
			double delta4E = phi_(ic     , jc + 2) - 3 * phi_(ic     , jc + 1) + 3 * phi_(ic      , jc    ) - phi_(ic     , jc - 1);
			double delta4N = phi_(ic + 1 , jc    ) - 3 * phi_(ic     , jc    ) + 3 * phi_(ic - 1  , jc    ) - phi_(ic - 2 , jc    );
			double delta4S = phi_(ic + 2 , jc    ) - 3 * phi_(ic + 1 , jc    ) + 3 * phi_(ic      , jc    ) - phi_(ic - 1 , jc    );

			// The signs should be corrected
			DCoeff_(i , j).w = s2_(ic , jc).w * xFaces(i , j).w * lambda_(ic , jc).w * delta2W - s4_(ic , jc).w * xFaces(i , j).w * lambda_(ic , jc).w * delta4W;
			
			DCoeff_(i , j).e = s2_(ic , jc).e * xFaces(i , j).e * lambda_(ic , jc).e * delta2E - s4_(ic , jc).e * xFaces(i , j).e * lambda_(ic , jc).e * delta4E;
			
			DCoeff_(i , j).n = s2_(ic , jc).n * yFaces(i , j).n * lambda_(ic , jc).n * delta2N - s4_(ic , jc).n * yFaces(i , j).n * lambda_(ic , jc).n * delta4N;

			DCoeff_(i , j).s = s2_(ic , jc).s * yFaces(i , j).s * lambda_(ic , jc).s * delta2E - s4_(ic , jc).s * yFaces(i , j).s * lambda_(ic , jc).s * delta4S;
			// DCoeff_(i , j).w =    s2_(ic , jc).w * std::abs(xFaces(i , j).w) * lambda_(ic , jc).w * (phi_(ic , jc) - phi_(ic , jc - 1)) -
			// 						    s4_(ic , jc).w * std::abs(xFaces(i , j).w) * lambda_(ic , jc).w * ( ( phi_(ic , jc + 2) - 2*phi_(ic , jc + 1) + phi_(ic , jc ) )- delta2Eta(phi_ , i , j - 1) ) ;
			// DCoeff_(i , j).e =    s2_(ic , jc).e * std::abs(xFaces(i , j).e) * lambda_(ic , jc).e * (phi_(ic , jc + 1) - phi_(ic , jc)) -
			// 					       s4_(ic , jc).e * std::abs(xFaces(i , j).e) * lambda_(ic , jc).e * ( delta2Eta(phi_ , i , j + 1) - delta2Eta(phi_ , i , j) );
			
			// DCoeff_(i , j).n =    s2_(ic , jc).n * std::abs(yFaces(i , j).n) * lambda_(ic , jc).n * (phi_(ic  , jc) - phi_(ic - 1, jc)) -
			// 					       s4_(ic , jc).n * std::abs(yFaces(i , j).n) * lambda_(ic , jc).n * ( delta2Csi(phi_ , i , j) - delta2Csi(phi_ , i - 1 , j) );
			// DCoeff_(i , j).s =    s2_(ic , jc).s * std::abs(yFaces(i , j).s) * lambda_(ic , jc).s * (phi_(ic + 1 , jc) - phi_(ic  , jc)) -
			// 						    s4_(ic , jc).s * std::abs(yFaces(i , j).s) * lambda_(ic , jc).s * ( delta2Csi(phi_ , i + 1, j) - delta2Csi(phi_ , i , j) );

			D_(i , j) =  DCoeff_(i , j).s - DCoeff_(i , j).n + DCoeff_(i , j).e - DCoeff_(i , j).w;
			// D_(i , j) =   DCoeff_(i , j).e - DCoeff_(i , j).w;
						
		}
	}
	
}

double Variable::computeError()
{
	double totalResidual(0.);
	double errMax = 0;
	for (unsigned int i = 0; i < Nci_; ++i)
	{
		for (unsigned int j = 1; j < Mci_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;

			totalResidual = totalResidual + (D_(i,j) - R_(i,j));
			errMax = max( std::abs(phi_(ic,jc) - phiPrev_(ic,jc)) , errMax)/phi_(ic, jc);

		}
	}

	return errMax;
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