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
	void sweepRows(char var);

	void buildColMatrix(Matrix<double>& Ay, unsigned int j);  
	void buildColVector_x(Vector<double>& by, unsigned int j);  
	void buildColVector_y(Vector<double>& bx, unsigned int j);  
	void sweepCols();

	Vector<double> GaussSeidel(Matrix<double>& A, Vector<double>& b, double tol=1e-5 , unsigned int N=10);
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
_deltaCsi(_L/(N+2)),
_deltaEta(_H/(M+2)),
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

	_y.setRow(_Jmax,  0 );
	// Adjust boundary conditions for y at i=0 and i = Imax
	_y.setRow( 0   , _H );
	for (int i = 1; i < _y.rows() - 1; ++i)
	{
		_y(i,0) = _y(i-1,0) - _deltaCsi;
		_y(i,_Jmax) = _y(i-1,_Jmax) - _deltaCsi;
	}


	_x.setColumn(_Imax, _L);
	// Adjust boundary conditions for x at j=0 and j= Jmax
	for (int j = 1; j < _x.cols() - 1; ++j)
	{
		_x(0,j) = _x(0,j-1) + _deltaEta;
		_x(_Imax,j) = _x(_Imax,j-1) + _deltaEta;
	}
	_x.setColumn(0  ,  0);

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

	std::cout << "y : " << std::endl;
	_y.print();
	std::cout << std::endl;


}

void Problem::solve()
{
	sweepRows('y');
	
	updateAlpha();
	updateBeta();
	updateGamma();
	updateCoeff();

	sweepCols();
	updateAlpha();
	updateBeta();
	updateGamma();
	updateCoeff();

	// sweepCols();
}


void Problem::sweepCols()
{
	std::cout << "Sweeping over the cols " << std::endl;

	// At each column of the domain (from j=1 to j=M_)
	// A matrix and a vector are created. Then the Gauss-Seidel algorithm is used.
	for (int j = _M; j > 0; --j)   // sweep over all jnternal rows. From j=1 to j=_M
	{
		Matrix<double> Axj(_N, _N); // Renamed Ay into Aj to highlight the dependency on the column
		Vector<double> bxj(_N);

		buildColMatrix(Axj, j);
		buildColVector_x(bxj, j);
		std::cout << "Axj : " << std::endl;
		Axj.print();

		std::cout << "bxj : " << std::endl;
		bxj.print();

		// Store the solution
		Vector<double> solxj = GaussSeidel(Axj,bxj);

		std::cout << "sol " << std::endl;
		solxj.print();

		// Insert solution inside the domain
		for (int i = 1; i < _N + 1; ++i)    // Insert solution
		{
			_x(i,j) = solxj(i-1);
		}


		
		std::cout << "x : " << std::endl;
		_x.print();

		updateAlpha();
		updateBeta();
		updateGamma();
		updateCoeff();

		

	}
}

void Problem::sweepRows(char var)
{
	// sweep over all internal rows. From i=1 to i=_N (i < _N + 1)
	// build a linear system in each row (matrix + vector)
	// solve the linear system and store the solution
	// Fill the domain with the solution
	if (var == 'x')
	{

		std::cout << "Sweeping over the rows " << std::endl;
		for (int i = 1; i < _N + 1; ++i)   
		{
			Matrix<double> Axi(_M, _M);   // These resources should be freed after the solution has been obtained
			Vector<double> bxi(_M);

			buildRowMatrix(Axi, i);
			buildRowVector_x(bxi, i);

			std::cout << "Ax : " << std::endl;
			Axi.print();

			std::cout << "bx : " << std::endl;
			bxi.print();

			Vector<double> solxi = GaussSeidel(Axi,bxi);
			std::cout << "solution : " <<  std::endl;
			solxi.print();

			for (int j = 1; j < _M + 1; ++j)    // Insert solution looping over internal points. The internal index j correspond to k = j - 1
			{
				_x(i,j) = solxi(j-1);
			}

			std::cout << "solx " << std::endl;
			solxi.print();
			std::cout << "x : " << std::endl;
			_x.print();
		}
	}

	if (var=='y')
	{
	
		for (int i = 1; i < _N + 1; ++i)
		{
			Matrix<double> Ayi(_M, _M);
			Vector<double> byi(_M);
			buildRowMatrix(Ayi, i);
			buildRowVector_y(byi, i);

			std::cout << "Ay : " << std::endl;
			Ayi.print();

			std::cout << "by : " << std::endl;
			byi.print();

			Vector<double> solyi = GaussSeidel(Ayi,byi);

			for (int j = 1; j < _M + 1; ++j)    // Insert solution
			{
				_y(i,j) = solyi(j-1);
			}

			std::cout << "soly " << std::endl;
			solyi.print();
			std::cout << "y : " << std::endl;
			_y.print();

			updateAlpha();
			updateBeta();
			updateGamma();
			updateCoeff();
		}

		_x.print();

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

	std::cout << "Updated gamma : " << std::endl;
	 _gamma.print() ;
}


void Problem::updateCoeff()
{
	for (int i = 0; i < _N + 2; ++i)
	{
		for (int j = 0; j < _M + 2; ++j)
		{
			_aP(i,j) = -2*_alpha(i,j)/pow(_deltaCsi, 2) - 2*_gamma(i,j)/pow(_deltaCsi, 2);
			_aN(i,j) = _alpha(i,j)/pow(_deltaCsi,2);
			_aS(i,j) = _alpha(i,j)/pow(_deltaCsi,2);
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

	std::cout << "aNE : " << std::endl;
	_aNE.print();
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

void Problem::buildRowVector_y(Vector<double>& by, unsigned int i)  // integer i locates the row of the domain
{
	for (unsigned int j = 1; j < by.size() + 1; ++j) // loop over internal points of the domain. Given a j, l is determined as l = j - 1. l is used to fill the b vector
	{
		unsigned int l = j - 1;
		by(l) = _aNE(i,j)*(_y(i+1, j+1) - _y(i+1, j-1) - _y(i-1, j+1) + _y(i-1, j-1) ) - _aN(i,j)*(_y(i+1,j) + _y(i-1,j)); 
	}
}



void Problem::buildColMatrix(Matrix<double>& A, unsigned int j)
{
	for (int l = 0; l < A.rows(); ++l)
	{
		int k = l + 1;
		A(l,l) = _aP(k,j);   // set diagonal. There should be an index that travels internally to the domain, like k = l +1
		if (l < _N - 1)
		{
			A(l,l+1) = _aS(k+1,j);
		}
		if (l > 0)
		{
			A(l,l-1) = _aN(k-1,j);
		}
	}

}

// temporary. I don't remember why.
void Problem::buildColVector_y(Vector<double>& by, unsigned int j)  // integer j locates the column of the domain
{
	for (unsigned int i = 1; i < by.size() + 1; ++i) // loop over internal points of the domain. Given a i, l is determined as l = i - 1. l is used to fill the b vector
	{
		unsigned int l = i - 1;
		std::cout << "_aE(i,j) = " << _aE(i,j) << std::endl;
		std::cout << "_x(i,j+1) = " << _x(i,j+1) << std::endl;
		std::cout << "_x(i,j-1) = " << _x(i,j-1) << std::endl;
		std::cout << "_aNE(i,j) = " << _aNE(i,j) << std::endl;
		std::cout << "x(i+1,j+1) = " << _x(i+1,j+1) << std::endl;
		std::cout << "x(i-1,j+1) = " << _x(i-1,j+1) << std::endl;
		std::cout << "x(i+1,j-1) = " << _x(i+1,j-1) << std::endl;
		std::cout << "x(i-1,j-1) = " << _x(i-1,j-1) << std::endl;
		std::cout << "_aNE(i,j)*(_x(i+1, j+1) - _x(i+1, j-1) - _x(i-1, j+1) + _x(i-1, j-1) ) = " << _aNE(i,j)*(_x(i+1, j+1) - _x(i+1, j-1) - _x(i-1, j+1) + _x(i-1, j-1) ) << std::endl;
		std::cout << "_aE(i,j)*(_x(i,j+1) + _x(i,j-1)) = " << _aE(i,j)*(_x(i,j+1) + _x(i,j-1)) << std::endl;
		by(l) = _aNE(i,j)*(_y(i+1, j+1) - _y(i+1, j-1) - _y(i-1, j+1) + _y(i-1, j-1) ) - _aE(i,j)*(_y(i,j+1) + _y(i,j-1)); 
		
	}
}

void Problem::buildColVector_x(Vector<double>& bx, unsigned int j)  // integer j locates the column of the domain
{
	std::cout << "x : " << std::endl;
	_x.print();
	for (unsigned int i = 1; i < bx.size() + 1; ++i) // loop over internal points of the domain. Given a i, l is determined as l = i - 1. l is used to fill the b vector
	{
		unsigned int l = i - 1;
		std::cout << "_aE(i,j) = " << _aE(i,j) << std::endl;
		std::cout << "_x(i,j+1) = " << _x(i,j+1) << std::endl;
		std::cout << "_x(i,j-1) = " << _x(i,j-1) << std::endl;
		std::cout << "_aNE(i,j) = " << _aNE(i,j) << std::endl;
		std::cout << "x(i+1,j+1) = " << _x(i+1,j+1) << std::endl;
		std::cout << "x(i-1,j+1) = " << _x(i-1,j+1) << std::endl;
		std::cout << "x(i+1,j-1) = " << _x(i+1,j-1) << std::endl;
		std::cout << "x(i-1,j-1) = " << _x(i-1,j-1) << std::endl;
		std::cout << "_aNE(i,j)*(_x(i+1, j+1) - _x(i+1, j-1) - _x(i-1, j+1) + _x(i-1, j-1) ) = " << _aNE(i,j)*(_x(i+1, j+1) - _x(i+1, j-1) - _x(i-1, j+1) + _x(i-1, j-1) ) << std::endl;
		std::cout << "_aE(i,j)*(_x(i,j+1) + _x(i,j-1)) = " << _aE(i,j)*(_x(i,j+1) + _x(i,j-1)) << std::endl;
		bx(l) = _aNE(i,j)*(_x(i+1, j+1) - _x(i+1, j+1) - _x(i+1, j+1) + _x(i-1, j-1) ) - _aE(i,j)*(_x(i,j+1) + _x(i,j-1)); 
		std::cout << "bx(l) = " << bx(l) << std::endl;
	}
}

Vector<double> Problem::GaussSeidel(Matrix<double>& A, Vector<double>& b, double tol , unsigned int N)
{
	Vector<double> u(A.cols());
    Vector<double> uPrev(A.cols());
    double err = tol + 1;  // Initialize err to a value greater than tol
    unsigned int n = 0;

    while (n < N && err > tol)
    {
        uPrev = u;  // Save the previous iteration

        for (int l = 0; l < A.cols(); ++l)
        {
            u(l) = b(l) / A(l, l);

            for (int k = 0; k < A.cols(); ++k)
            {
                if (l != k)
                {
                    u(l) -= A(l, k) / A(l, l) * u(k);
                }
            }
        }

        // Calculate the error
        err = 0.0;
        for (int i = 0; i < A.cols(); ++i)
        {
            err += std::abs(u(i) - uPrev(i));
        }

        ++n;

        if (n > N)
        {
        	std::cout << "CONVERGENCE FAILED " << std::endl;
        }
    }

    return u;
	

	return u;
}
