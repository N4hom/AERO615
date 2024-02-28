#include <iostream>
#include <math.h>
#include <cmath>
#include "Matrix.hpp"
#include "Vector.hpp"

class Problem
{
private:

	double _L = 1.0, _H = 1.0;   // Lenght and height
	unsigned int _N, _M;  // Size of the problem
	unsigned int _Imax;
	unsigned int _Jmax;
	double _deltaCsi, _deltaEta;
	Matrix<double> _x, _y;  // Solution storage
	Matrix<double> _alpha, _beta, _gamma;  // Coefficient storage. They could be a vector but for now I'll keep them Matrices. I was kidding, they're matrices
	Matrix<double> _aP, _aE, _aW, _aN, _aS, _aNE, _aSE, _aSW, _aNW;   // Storage for linear system coefficient

	Matrix<double> _A; // Problem matrix


public:
	Problem(unsigned int, unsigned int);
	~Problem();
	void initialize();
	void updateCoeff();
	void updateAlpha();
	void updateBeta();
	void updateGamma();
	
	// Create linear system in each row of the domain. To be used inside sweepRows() where a temporary matrix is created
	void buildRowMatrix(Matrix<double>& Ax, unsigned int i);  
	void buildRowVector_x(Vector<double>& bx, unsigned int i);  
	void buildRowVector_y(Vector<double>& by, unsigned int i);  
	void sweepRows();

	void buildColMatrix(Matrix<double>& Ay, unsigned int j);  
	void buildColVector(Vector<double>& by, unsigned int j);  
	void sweepCols();

	Vector<double> GaussSeidel(Matrix<double>& A, Vector<double>& b, double tol=1e-5 , unsigned int N=100);
	void solve();

	Matrix<double>& alpha(){ return _alpha;};
	Matrix<double>& beta(){ return _beta;};
	Matrix<double>& gamma(){ return _gamma;};

	Matrix<double>& x(){ return _x;};
	Matrix<double>& y(){ return _y;};

	
};

Problem::Problem(unsigned int N, unsigned int M):
_N(N),
_M(M),
_Imax(N + 1),
_Jmax(M + 1),
_deltaCsi(_L/N),
_deltaEta(_H/M),
_x(N + 2 , M + 2),
_y(N + 2 , M + 2),
_alpha(N + 2 , M + 2),
_beta(N + 2 , M + 2),
_gamma(N + 2 , M + 2),
_aP(N + 2 , M + 2),
_aE(N + 2 , M + 2),
_aW(N + 2 , M + 2),
_aN(N + 2 , M + 2),
_aS(N + 2 , M + 2),
_aNE(N + 2 , M + 2),
_aSE(N + 2 , M + 2),
_aSW(N + 2 , M + 2),
_aNW(N + 2 , M + 2),
_A(N + 2 , M + 2)
{
	std::cout << "Problem initialized " << std::endl;
	std::cout << "_deltaEta =" << _deltaEta << std::endl;
	std::cout << "_deltaCsi " << _deltaCsi << std::endl;
	std::cout << "Number of nodes along x : " << N << std::endl;
	std::cout << "Number of nodes along y : " << M << std::endl;
}

Problem::~Problem()
{
	std::cout << "Problem over " << std::endl;
}

void Problem::initialize()
{
	// x and y are set to zero everywhere except on the boundary 

	_y.setRow( 0   , _H );
	_y.setRow(_Jmax,  0 );

	_x.setColumn(0  ,  0);
	_x.setColumn(_Imax, _L);

	// alpha, beta and gamma are set to zero everywhere except in the interior

	updateAlpha();
	updateBeta();
	updateGamma();

	// Initialize the a coefficients. 9 Matrices to be created.
	updateCoeff();

	std::cout << "alpha : " << std::endl;
	_alpha.print();
	std::cout << std::endl;

	std::cout << "beta : " << std::endl;
	_beta.print();
	std::cout << std::endl;

	std::cout << "gamma : " << std::endl;
	_gamma.print();
	std::cout << std::endl;


	std::cout << "x : " << std::endl;
	_x.print();
	std::cout << std::endl;


}

void Problem::solve()
{
	sweepRows();
	
	updateAlpha();
	updateBeta();
	updateGamma();
	updateCoeff();

	// sweepCols();
}


void Problem::sweepCols()
{
	for (int j = _M; j > 0; --j)   // sweep over all jnternal rows. From j=1 to j=_M
	{
		Matrix<double> Ay(_N, _N); // It should be named Ax because the x axis is vertical in my ref
		Vector<double> by(_N);

		buildColMatrix(Ay, j);
		buildColVector(by, j);
		_x.setColumn(j, GaussSeidel(Ay,by));

		std::cout << "Ay : " << std::endl;
		Ay.print();

		std::cout << "by : " << std::endl;
		by.print();

		
		std::cout << "x : " << std::endl;
		_x.print();

		// Ay.print();
		// std::cout << "b: " << std::endl;
		// by.print();
		// std::cout << std::endl;

	}
}

void Problem::sweepRows()
{
	for (int i = 1; i < _N + 1; ++i)   // sweep over all internal rows. From i=1 to i=_N
	{
		Matrix<double> Ax(_M, _M);
		Vector<double> bx(_M);

		buildRowMatrix(Ax, i);
		buildRowVector_x(bx, i);

		std::cout << "Ax : " << std::endl;
		Ax.print();

		std::cout << "bx : " << std::endl;
		bx.print();

		Vector<double> solx = GaussSeidel(Ax,bx);
		std::cout << "solution : " <<  std::endl;
		solx.print();


		Matrix<double> Ay(_M, _M);
		Vector<double> by(_M);
		buildRowMatrix(Ay, i);
		buildRowVector_y(by, i);

		std::cout << "Ay : " << std::endl;
		Ax.print();

		std::cout << "by : " << std::endl;
		bx.print();

		Vector<double> soly = GaussSeidel(Ay,by);

		for (int j = 1; j < _M; ++j)    // Insert solution
		{
			_y(i,j) = soly(j);
		}

		std::cout << "soly " << std::endl;
		soly.print();
		std::cout << "y : " << std::endl;
		_y.print();

	}
}

void Problem::updateAlpha()
{
	for (int i = 1; i < _alpha.rows() - 1; ++i)
	{
		for (int j = 1; j < _alpha.cols() - 1; ++j)
		{
			
			_alpha(i,j) = 1.0/4.0/pow(_deltaEta, 2) *(pow( _x(i,j + 1) - _x(i,j - 1) , 2) + pow(_y(i,j + 1) - _y(i,j-1) , 2));
			
		}
	}
}

void Problem::updateBeta()
{

	for (int i = 1; i < _alpha.rows() - 1; ++i)
	{
		for (int j = 1; j < _alpha.cols() - 1; ++j)
		{
			_beta(i,j) = 1.0/4./_deltaEta/_deltaCsi*((_x(i+1,j) - _x(i-1,j))*(_x(i,j+1) - _x(i, j -1)) + (_y(i+1,j - _y(i-1,j)))*(_y(i,j+1) - _y(i,j-1)) );	
		}
	}
}


void Problem::updateGamma()
{
	for (int i = 1; i < _alpha.rows() - 1; ++i)
	{
		for (int j = 1; j < _alpha.cols() - 1; ++j)
		{
			_gamma(i,j) = 1.0/4.0/pow(_deltaCsi, 2) *(pow( _x(i + 1,j ) - _x(i - 1,j ) , 2) + pow(_y(i + 1,j) - _y(i - 1,j) , 2));	
		}
	}
}


void Problem::updateCoeff()
{
	for (int i = 0; i < _N + 2; ++i)
	{
		for (int j = 0; j < _M + 2; ++j)
		{
			_aP(i,j) = -2*_alpha(i,j)/pow(_deltaCsi, 2) - 2*_gamma(i,j)/pow(_deltaCsi, 2);
			_aN(i,j) = 2*_alpha(i,j)/pow(_deltaCsi,2);
			_aS(i,j) = 2*_alpha(i,j)/pow(_deltaCsi,2);
			_aNE(i,j) = _beta(i,j)/2/_deltaCsi/_deltaEta;
			_aSE(i,j) = _beta(i,j)/2/_deltaCsi/_deltaEta;
			_aSW(i,j) = _beta(i,j)/2/_deltaCsi/_deltaEta;
			_aNW(i,j) = _beta(i,j)/2/_deltaCsi/_deltaEta;
			_aE(i,j) = _gamma(i,j)/pow(_deltaEta,2);
			_aW(i,j) = _gamma(i,j)/pow(_deltaEta,2);
		}
	}

	std::cout << "aE : " << std::endl;
	_aE.print();

	std::cout << "aP : " << std::endl;
	_aP.print();

	std::cout << "aW : " << std::endl;
	_aW.print();
}


void Problem::buildRowMatrix(Matrix<double>& A, unsigned int i)
{
	for (int l = 0; l < A.cols(); ++l)
	{
		int k = l + 1;
		A(l,l) = _aP(i,k);   // set diagonal. There should be an index that travels internally to the domain, like k = l +1
		if (l < _N - 1)
		{
			A(l,l+1) = _aE(i,k+1);
		}
		if (l > 0)
		{
			A(l,l-1) = _aW(i,k-1);
		}
	}

}

void Problem::buildRowVector_x(Vector<double>& bx, unsigned int i)  // integer i locates the row of the domain
{
	for (unsigned int j = 1; j < bx.size() + 1; ++j) // loop over internal points of the domain. Given a j, l is determined as l = j - 1. l is used to fill the b vector
	{
		unsigned int l = j - 1;
		bx(l) = _aNE(i,j)*(_x(i+1, j+1) - _x(i+1, j-1) - _x(i-1, j+1) + _x(i-1, j-1) ) - _aN(i,j)*(_x(i+1,j) + _x(i-1,j)); 
	}
}

void Problem::buildRowVector_y(Vector<double>& bx, unsigned int i)  // integer i locates the row of the domain
{
	for (unsigned int j = 1; j < bx.size() + 1; ++j) // loop over internal points of the domain. Given a j, l is determined as l = j - 1. l is used to fill the b vector
	{
		unsigned int l = j - 1;
		bx(l) = _aNE(i,j)*(_y(i+1, j+1) - _y(i+1, j-1) - _y(i-1, j+1) + _y(i-1, j-1) ) - _aN(i,j)*(_y(i+1,j) + _y(i-1,j)); 
	}
}



void Problem::buildColMatrix(Matrix<double>& A, unsigned int j)
{
	for (int l = 0; l < A.rows(); ++l)
	{
		A(l,l) = _aP(l,j);   // set diagonal. It should be _aP(l+1,j)
		if (l < _N - 1)
		{
			A(l,l+1) = _aS(l+1,j);
		}
		if (l > 0)
		{
			A(l,l-1) = _aN(l-1,j);
		}
	}

}

// temporary
void Problem::buildColVector(Vector<double>& by, unsigned int j)  // integer j locates the column of the domain
{
	for (unsigned int i = 1; i < by.size() + 1; ++i) // loop over internal points of the domain. Given a i, l is determined as l = i - 1. l is used to fill the b vector
	{
		unsigned int l = i - 1;
		by(l) = _aNE(i,j)*(_x(i+1, j+1) - _x(i+1, j-1) - _x(i-1, j+1) + _x(i-1, j-1) ) - _aE(i,j)*(_x(i,j+1) + _x(i,j-1)); 
		if (j==_M)
		{
			std::cout << "by(l) for l = " << l << " " << by(l) << std::endl;
			std::cout << "_x(i,j+1) = " << _x(i,j+1) <<  std::endl;
			std::cout << "_x(i,j-1) = " << _x(i,j-1) <<  std::endl;
		}
	}
}

Vector<double> Problem::GaussSeidel(Matrix<double>& A, Vector<double>& b, double tol , unsigned int N)
{
	Vector<double> u = b;
	Vector<double> uPrev(A.cols());
	double err = 1;
	unsigned int n = 1;

	while(n < N)
	{
		for (int l = 0; l < A.cols() ; ++l)
		{
			u(l) = b(l)/A(l,l);

			for (int k = 0; k < A.cols(); ++k)
			{
				if (l != k)
				{
					u(l) -= A(l,k)/A(l,l) * u(k);
				}
			}

			if (std::isnan(u(l))) 
			{
			    u(l) = 0.;
			}
		}

		++n;

	}

	

	return u;
}
