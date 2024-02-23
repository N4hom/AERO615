#include <iostream>
#include <math.h>
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
	Matrix<double> _A; // Problem matrix
public:
	Problem(unsigned int, unsigned int);
	~Problem();
	void initialize();
	void updateCoeff();
	void updateAlpha();
	void updateBeta();
	void updateGamma();
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
_A(N + 2 , M + 2)
{
	std::cout << "Problem initialized " << std::endl;
	std::cout << "_deltaEta " << _deltaEta << std::endl;
	std::cout << "pow(_deltaEta, 2) " << pow(_deltaEta, 2) << std::endl;
	std::cout << "_deltaCsi " << _deltaCsi << std::endl;
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

	// update alpha
	// update beta 
	// update gamma
	// update coeff

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
			//std::cout << _x(i, j + 1 ) - _x(i, j - 1 )<< std::endl;
			//std::cout << pow( _x(i,j + 1) - _x(i,j - 1) , 2) + pow(_y(i,j + 1) - _y(i,j-1) , 2)<< std::endl;
			std::cout << (pow( _x(i + 1,j ) - _x(i - 1,j ) , 2) + pow(_y(i + 1,j) - _y(i - 1,j) , 2)) << std::endl;
			_gamma(i,j) = 1.0/4.0/pow(_deltaCsi, 2) *(pow( _x(i + 1,j ) - _x(i - 1,j ) , 2) + pow(_y(i + 1,j) - _y(i - 1,j) , 2));
			
		}
	}
}
