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

	// Adjust boundary conditions for y at i=0 and i = Imax

	std::cout << "Initializing x " << std::endl;

	_x.setColumn(_Jmax, _L);
	// Adjust boundary conditions for x at j=0 and j= Jmax
	for (int j = 1; j < _x.cols() - 1; ++j)
	{
		_x(0,j) = _x(0,j-1) + _deltaEta;
		_x(_Imax,j) = _x(_Imax,j-1) + _deltaEta;

		for (int i = 1; i < _x.rows() - 1; ++i)
		{
			//_x(i,j) = 0.0 ;
			_x(i,j) = _x(i-1,j) ;
		}
	}
	_x.setColumn(0  ,  0);

	std::cout << "Initializing y " << std::endl;
	

	for (int j = 0; j < _y.cols() ; ++j)
	{
		double etaj = j * _deltaEta;

		if (etaj > 2 && etaj < 3)
		{
			_y(0,j) = 1 - 0.2*sin(M_PI *(_x(0,j) - 2));
			_y(_Imax,j) = 0.2*sin(M_PI *(_x(0,j) - 2));

		}
		else
		{
			_y(0,j) = _H;
			_y(_Imax,j) = 0;
		}
	}



	for (int i = 1; i < _y.rows() - 1; ++i)
	{
		_y(i,0) = _y(i-1,0) - _deltaCsi;
		_y(i,_Jmax) = _y(i-1,_Jmax) - _deltaCsi;

	}

	std::cout << "Initial x " << std::endl;
	_x.print();
	std::cout << "Initial y " << std::endl;
	_y.print();

}

void Problem::solve2()
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

		
		for (int i = 1; i < _N + 1; ++i)
		{
			for (int j = 1; j < _M + 1; ++j)
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

		for (int j = _M ; j > 0 ; --j)
		{
			for (int i = 1; i < _N + 1; ++i)
			{
				updateAlpha(i,j);
				updateBeta(i,j);
				updateGamma(i,j);
				_a1 =  _beta/2/_deltaCsi/_deltaEta/-2.0/(_alpha/pow2(_deltaCsi) + _gamma/pow2(_deltaEta));
				_a2 = -_gamma/pow2(_deltaEta)/-2.0/(_alpha/pow2(_deltaCsi) + _gamma/pow2(_deltaEta));
				_a4 = -_alpha/pow2(_deltaCsi)/-2.0/(_alpha/pow2(_deltaCsi) + _gamma/pow2(_deltaEta));

				_x(i,j) = _a1 * ( _x(i+1,j+1) - _x(i-1,j+1) - _x(i+1,j-1) + _x(i-1, j-1) ) + _a2 * (_x(i,j+1) + _x(i,j-1)) + _a4 * (_x(i+1,j) + _x(i-1,j));


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

		if (iter > 1000)
		{
			std::cout << "Not converged !!!!!!!!!!!!!" << std::endl;
			std::cout << "Not converged !!!!!!!!!!!!!" << std::endl;
			std::cout << "Not converged !!!!!!!!!!!!!" << std::endl;
			std::cout << "Not converged !!!!!!!!!!!!!" << std::endl;
			break;
		}
		// _x.print();
		// _y.print();

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


// void Problem::updateCoeff()
// {
// 	std::cout << "Updating coeff " << std::endl;
// 	_a1 = _beta(i,j)/2/_deltaCsi/_deltaEta/-2.0/(_alpha(i,j)/_deltaCsi + _gamma(i,j)/_deltaEta);
// 	_a2 = - _gamma(i,j)/pow2(_deltaEta)/-2/(_alpha(i,j)/_deltaCsi + _gamma(i,j)/_deltaEta);
// 	_a4 = - _alpha(i,j)/pow2(_deltaEta)/-2/(_alpha(i,j)/_deltaCsi + _gamma(i,j)/_deltaEta);
	
// }




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
        	std::cout << "CONVERGENCE FAILED !!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        }
    }

    return u;
	

	return u;
}

// double calculateError()
// {
// 	for (int i = 1; i < _N - 1 ; ++i)
// 	{
// 		for (int j = 1; j < _M + 1; ++i)
// 		{
// 			_yPrev
// 		}
// 	}
// }


void Problem::save()
{
	_y.save("y");
	_x.save("x");
}

