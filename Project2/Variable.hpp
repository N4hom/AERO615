#include "Matrix.hpp"
#include "Vector.hpp"

class Variable
{

	unsigned int    Nci_;
	unsigned int    Mci_;
	unsigned int    Nc_;
	unsigned int    Mc_;
	Matrix<double>  phi_;
	Matrix<double>  flux_f_;
	Matrix<double>  flux_g_;
	Matrix<double>  R_;
	Matrix<double>  D_;
	Mesh& 			mesh;

public:
	Variable(unsigned int N , unsigned int M , Mesh& mesh);
	~Variable();


	double interpolateLeft  (Matrix<double>& flux , unsigned int i, unsigned int j) const;
	double interpolateRight (Matrix<double>& flux , unsigned int i, unsigned int j) const;
	double interpolateTop   (Matrix<double>& flux , unsigned int i, unsigned int j) const;
	double interpolateBottom(Matrix<double>& flux , unsigned int i, unsigned int j) const;

	double flux(Matrix<double>& velocity , unsigned int i, unsigned int j);

	void computeFlux_f(Matrix<double>& velocity                                 );
	void computeFlux_g(Matrix<double>& velocity                                 );
	void correctFlux_f(Matrix<double>& velocity , unsigned int i, unsigned int j);
	void correctFlux_g(Matrix<double>& velocity , unsigned int i, unsigned int j);

	double computeResidualij(unsigned int i , unsigned int j                    );
	void computeResidual();
	double Rij(unsigned int i, unsigned int j) const{return R_(i,j);};
	Matrix<double>& R(){return R_;};
	Matrix<double>& D(){return D_;};


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
Variable::Variable(unsigned int N, unsigned int M , Mesh& mesh):
Nci_(N - 4),
Mci_(M - 4),
Nc_(N),
Mc_(M),
phi_(Nc_ , Mc_),
flux_f_(Nc_ , Mc_),
flux_g_(Nc_ , Mc_),
R_(Nci_ , Mci_),
D_(Nci_ , Mci_),
mesh(mesh)
{
	std::cout << "R " << std::endl;
	R_.print();
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
			flux_f_(i , j) = flux(velocity , ic , jc);
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
			flux_g_(i , j) = flux(velocity , ic , jc);
		}
	}
}

void Variable::correctFlux_g(Matrix<double>& velocity , unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;

	flux_g_(i , j) = flux(velocity , ic , jc);

}

void Variable::correctFlux_f(Matrix<double>& velocity , unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	
	flux_f_(i , j) = flux(velocity , ic , jc);

}


double Variable::computeResidualij(unsigned int i, unsigned int j)
{
	unsigned int ic = i + 2;
	unsigned int jc = j + 2;
	
	double RfLeft = interpolateLeft(flux_f_ , ic , jc) * mesh.yFacesLeft_(i , j);  // f component calculated assuming that yFacesLeft_ stores face area along csi(y)
	double RfRight = interpolateRight(flux_f_ , ic , jc) * mesh.yFacesRight_(i , j);
	double RfTop = interpolateTop(flux_f_ , ic , jc) * mesh.yFacesTop_(i , j);;
	double RfBottom = interpolateBottom(flux_f_ , ic , jc) * mesh.yFacesBottom_(i , j);;

	double RgLeft = interpolateLeft(flux_g_ , ic , jc) * mesh.yFacesLeft_(i , j);;
	double RgRight = interpolateRight(flux_g_ , ic , jc) * mesh.yFacesRight_(i , j);;
	double RgTop = interpolateTop(flux_g_ , ic , jc) * mesh.yFacesTop_(i , j);;
	double RgBottom = interpolateBottom(flux_g_ , ic , jc) * mesh.yFacesBottom_(i , j);;


	return (RfLeft + RfRight + RfTop + RfBottom) - (RgLeft + RgRight + RgTop + RgBottom);
	
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

