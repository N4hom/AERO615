#include "Problem.hpp"
#include <math.h>
#include <cmath>


Problem::Problem(unsigned int N, unsigned int M):
_N(N),
_M(M),
_Imax(N + 1),
_Jmax(M + 1),
_deltaCsi(_H/(N+1)),
_deltaEta(_L/(M+1)),
_x(N + 2 , M + 2),
_y(N + 2 , M + 2),
// _xPrev(_N+2,_M+2),
// _yPrev(_N+2,_M+2),
_alpha(0.0),
_beta(0.0),
_gamma(0.0)
{
	std::cout << "Problem parameters " << std::endl;
	std::cout << "_deltaEta = " << _deltaEta << std::endl;
	std::cout << "_deltaCsi = " << _deltaCsi << std::endl;
	std::cout << "Number of nodes along csi : " << N << std::endl;
	std::cout << "Number of nodes along eta : " << M << std::endl;
}

Problem::~Problem()
{
}

double pow2(double f)
{
	return f*f;
}

void Problem::initialize()
{
	// x and y are set to zero everywhere except on the boundary 
	// Adjust boundary conditions for j = _Jmax ---> x = _L
	// and for j = 0 ---> x = 0.
	// Initialize internal points with algebraic grid

	std::cout << "Initializing x " << std::endl;

	_x.setColumn(_Jmax, _L);
	bool leadingEdge = true;
	bool trailingEdge = true;
		
	for (unsigned int j = 1; j < _x.cols() - 1; ++j)
	{

		double etaj = j * _deltaEta;

		_x(0,j) = _x(0,j-1) + _deltaEta;
		_x(_Imax,j) = _x(_Imax,j-1) + _deltaEta;
		
		if (etaj > 2 && leadingEdge)
		{

			_x(0,j - 1) = 2;
			_x(_Imax,j - 1) = 2;
			leadingEdge = false;
		}

		if (etaj > 3 && trailingEdge)
		{
			_x(0,j - 1) = 3;
			_x(_Imax,j - 1) = 3;
			trailingEdge = false;
		}
		

		for (unsigned int i = 1; i < _x.rows() - 1; ++i)
		{
			//_x(i,j) = 0.0 ;
			_x(i,j) = _x(i-1,j) ;
		}
	}
	_x.setColumn(0  ,  0);

	std::cout << "Initializing y " << std::endl;
	
	//- Applying boundary conditions at i = 0 and i = _Imax 
	//  if eta location is between  2 and 3 apply sinusoidal functions
	//  else set _y(0,j) to 0 and _y(_Imax,j) to _H

	for (unsigned int j = 0; j < _y.cols() ; ++j)
	{
		double etaj = j * _deltaEta;

		if (etaj > 2 && etaj < 3)
		{
			_y(0,j) = 1 - 0.1*sin(M_PI *(_x(0,j) - 2));  // f' = -0.2*M_PI*cos(M_PI *(_x(0,j) - 2))
			_y(_Imax,j) = 0.1*sin(M_PI *(_x(0,j) - 2));  // f' = 0.2*M_PI*cos(M_PI *(_x(0,j) - 2))

			//_y(0,j) = 1 ; 
			//_y(_Imax,j) = 0;

		}
		else
		{
			_y(0,j) = _H;
			_y(_Imax,j) = 0;
		}
	}


	//- Apply boundary conditions at i = 0  and i = _Imax for y
	for (unsigned int i = 1; i < _y.rows() - 1; ++i)
	{
		_y(i,0) = _y(i-1,0) - _deltaCsi;
		_y(i,_Jmax) = _y(i-1,_Jmax) - _deltaCsi;

	}

	std::cout << "Initial x " << std::endl;
	_x.print();
	std::cout << "Initial y " << std::endl;
	_y.print();

}

void Problem::solve()
{
	unsigned int iter = 0;
	Matrix<double> xPrev(_N+2,_M+2, 0);
	Matrix<double> yPrev(_N+2,_M+2, 0);
	double errX;
	double errY;

	do
	{	

		double errMax_x = 0;
		double errMax_y =  0;

		//- Sweep over rows for y
		for (unsigned int i = 1; i < _N + 1; ++i)
		{
			for (unsigned int j = 1; j < _M + 1; ++j)
			{

				updateAlpha(i,j);
				updateBeta(i,j);
				updateGamma(i,j);
				_a1 = _beta/2/_deltaCsi/_deltaEta/-2.0/(_alpha/pow2(_deltaCsi) + _gamma/pow2(_deltaEta));
				_a2 = - _gamma/pow2(_deltaEta)/-2.0/(_alpha/pow2(_deltaCsi) + _gamma/pow2(_deltaEta));
				_a4 = - _alpha/pow2(_deltaCsi)/-2.0/(_alpha/pow2(_deltaCsi) + _gamma/pow2(_deltaEta));

				_y(i,j) = _a1 * ( _y(i+1,j+1) - _y(i-1,j+1) - _y(i+1,j-1) + _y(i-1, j-1) ) + _a2 * (_y(i,j+1) + _y(i,j-1)) + _a4 * (_y(i+1,j) + _y(i-1,j));

				
				errMax_y = max( std::abs(_y(i,j) - yPrev(i,j)) , errMax_y);
				yPrev(i,j) = _y(i,j);

			}
		}

		//- Sweep over columns for x
		for (unsigned int j = _M ; j > 0 ; --j)
		{
			for (unsigned int i = 1; i < _N + 1; ++i)
			{
				updateAlpha(i,j);
				updateBeta(i,j);
				updateGamma(i,j);
				_a1 =  _beta/2/_deltaCsi/_deltaEta/-2.0/(_alpha/pow2(_deltaCsi) + _gamma/pow2(_deltaEta));
				_a2 = -_gamma/pow2(_deltaEta)/-2.0/(_alpha/pow2(_deltaCsi) + _gamma/pow2(_deltaEta));
				_a4 = -_alpha/pow2(_deltaCsi)/-2.0/(_alpha/pow2(_deltaCsi) + _gamma/pow2(_deltaEta));

				_x(i,j) = _a1 * ( _x(i+1,j+1) - _x(i-1,j+1) - _x(i+1,j-1) + _x(i-1, j-1) ) + _a2 * (_x(i,j+1) + _x(i,j-1)) + _a4 * (_x(i+1,j) + _x(i-1,j));

				double etaj = j * _deltaEta;

				//double fprime_0 = -0.2*M_PI*sin(M_PI *(_x(0,j) - 2));
				//double fprime_Imax = 0.2*M_PI*sin(M_PI *(_x(0,j) - 2));

				//Adjust boundary conditions for x
				if (etaj > 2 && etaj < 3)
				{
					// Top
					_x(1,j) = _x(0,j) - (-0.1*M_PI*cos(M_PI *(_x(0,j) - 2)))*(_y(1,j) - _y(0,j)) ;

					// Bottom
					_x(_Imax-1,j) = _x(_Imax,j) + (0.1*M_PI*cos(M_PI *(_x(0,j) - 2)))*(_y(_Imax,j) - _y(_Imax-1,j)) ;

				}
				

				errMax_x = max( std::abs(_x(i,j) - xPrev(i,j)) , errMax_x);
				xPrev(i,j) = _x(i,j);
			}
		}

		errX = errMax_x;
		errY = errMax_y;


		std::cout << "-------------------------------------" << std::endl;
		std::cout << "iter = " << iter << std::endl;
		std::cout << "max error for x : " << errX << std::endl;
		std::cout << "max error for y : " << errY << std::endl;
		std::cout << "-------------------------------------" << std::endl;

		++iter;

		if (iter > 100000)
		{
			std::cout << "Not converged !!!!!!!!!!!!!!!!!!!!!!" << std::endl;
			break;
		}
		

	}	
	while (errX > _tol || errY > _tol);


	
	
}





void Problem::updateAlpha(unsigned int i, unsigned int j)
{
	//std::cout << "Updating alpha " << std::endl;
	
			
	_alpha = 1.0/4.0/pow2(_deltaEta) *( pow2( _x(i,j + 1) - _x(i,j - 1)) + pow2(_y(i,j + 1) - _y(i,j-1)) );
			
	//std::cout << "alpha updated" << std::endl;

}

void Problem::updateBeta(unsigned int i, unsigned int j)
{
	//std::cout << "Updating beta " << std::endl;
	
			
	_beta = 1.0/4./_deltaEta/_deltaCsi*( (_x(i+1,j) - _x(i-1,j))*(_x(i,j+1) - _x(i, j -1)) + (_y(i+1,j) - _y(i-1,j))*(_y(i,j+1) - _y(i,j-1)) );	
		
	
	//std::cout << "Updated beta : " << std::endl;

}


void Problem::updateGamma(unsigned int i, unsigned int j)
{
	//std::cout << "Updating gamma " << std::endl;
	
	_gamma = 1.0/4.0/pow2(_deltaCsi) *(pow2( _x(i + 1,j ) - _x(i - 1,j )) + pow2(_y(i + 1,j) - _y(i - 1,j)));	
	
	//std::cout << "Updated gamma : " << std::endl;
}



void Problem::save()
{
	_y.save("y");
	_x.save("x");


	writeVTK(_x, _y, "output.vtk");
}

