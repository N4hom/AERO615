#include "Matrix.hpp"
#include "Vector.hpp"

class Variable
{

	unsigned int Nci_;
	unsigned int Mci_;
	unsigned int Nc_;
	unsigned int Mc_;
	Matrix<double> phi_;
	Matrix<double> flux_f_;
	Matrix<double> flux_g_;

public:
	Variable(unsigned int N , unsigned int M);
	~Variable();


	double interpolateLeft  (Matrix<double>& flux , unsigned int i, unsigned int j) const;
	double interpolateRight (Matrix<double>& flux , unsigned int i, unsigned int j) const;
	double interpolateTop   (Matrix<double>& flux , unsigned int i, unsigned int j) const;
	double interpolateBottom(Matrix<double>& flux , unsigned int i, unsigned int j) const;

	double flux(Matrix<double>& velocity , unsigned int i, unsigned int j);

	void computeFlux_f(Matrix<double>& velocity);
	void computeFlux_g(Matrix<double>& velocity);


	// Const access to the ith,jth element  (read only)
	const double &operator()(unsigned int i, unsigned int j) const;

	// Non const access for getting the ith, jth element
	double &operator()(unsigned int i, unsigned int j);

	void print() const{ phi_.print();};

	Matrix<double>& phi(){return phi_;}
	Matrix<double>& flux_f(){return flux_f_;}
	Matrix<double>& flux_g(){return flux_g_;}
	
};

// input is the number of total cell. Easier to set the Variable input using the size of the matrix phi
Variable::Variable(unsigned int N, unsigned int M):
Nci_(N - 4),
Mci_(M - 4),
Nc_(N),
Mc_(M),
phi_(Nc_ , Mc_),
flux_f_(Nci_ , Mci_),
flux_g_(Nci_ , Mci_)
{

}

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
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	return 0.5 * (flux(ic     , jc - 1) + flux(ic , jc)) ;
}

double Variable::interpolateRight(Matrix<double>& flux , unsigned int i, unsigned int j) const
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	return 0.5 * (flux(ic     , jc + 1) + flux(ic , jc)) ;
}

double Variable::interpolateTop(Matrix<double>& flux , unsigned int i, unsigned int j)	const
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	return 0.5 * (flux(ic + 1  , jc   ) + flux(ic , jc)) ;
}

double Variable::interpolateBottom(Matrix<double>& flux , unsigned int i, unsigned int j) const
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	return 0.5 * (flux(ic - 1  , jc   ) + flux(ic , jc)) ;
}

double Variable::flux(Matrix<double>& velocity , unsigned int i , unsigned int j)
{
	return phi_(i , j) * velocity(i , j);	
}

void Variable::computeFlux_f(Matrix<double>& velocity)
{
	for (unsigned int i = 0; i < Nc_; ++i)
	{
		for (unsigned int j = 0; j < Mc_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;
			// The fluxes are stored only in internal cells
			flux_f_(i , j) = phi_(ic , jc) * velocity(ic ,jc);
		}
	}
}

void Variable::computeFlux_g(Matrix<double>& velocity)
{
	for (unsigned int i = 0; i < Nc_; ++i)
	{
		for (unsigned int j = 0; j < Mc_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;
			// The fluxes are stored only in internal cells
			flux_g_(i , j) = phi_(ic , jc) * velocity(ic ,jc);
		}
	}
}

