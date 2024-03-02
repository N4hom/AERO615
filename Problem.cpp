#include "Problem.hpp"



Problem::Problem(unsigned int N, unsigned int M):
_N(N),
_M(M),
_Imax(N + 1),
_Jmax(M + 1),
_deltaCsi(_L/(N+1)),
_deltaEta(_H/(M+1)),
_x(N + 2 , M + 2),
_y(N + 2 , M + 2),
_alpha(N + 2 , M + 2),
_beta(N + 2 , M + 2),
_gamma(N + 2 , M + 2),
_aP(0.),
_aEW(0.),
_aNS(0.),
_aC(0.),
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

double pow2(double f)
{
	return f*f;
}

void Problem::initialize()
{
	// x and y are set to zero everywhere except on the boundary 

	// Adjust boundary conditions for y at i=0 and i = Imax
	_y.setRow(_Jmax,  0 );
	_y.setRow( 0   , _H );
	for (int i = 1; i < _y.rows() - 1; ++i)
	{

		_y(i,0) = _y(i-1,0) - _deltaCsi;
		_y(i,_Jmax) = _y(i-1,_Jmax) - _deltaCsi;

		for (int j = 1; j < _y.cols() - 1; ++j)
		{
//			_y(i,j) = _y(i,j-1) ;
			_y(i,j) = 0.0 ;
		}
	}


	_x.setColumn(_Imax, _L);
	// Adjust boundary conditions for x at j=0 and j= Jmax
	for (int j = 1; j < _x.cols() - 1; ++j)
	{
		_x(0,j) = _x(0,j-1) + _deltaEta;
		_x(_Imax,j) = _x(_Imax,j-1) + _deltaEta;

		for (int i = 1; i < _x.rows() - 1; ++i)
		{
			//_x(i,j) = _x(i-1,j) ;
			_x(i,j) = 0.0 ;
		}
	}
	_x.setColumn(0  ,  0);

	// alpha, beta and gamma are set to zero everywhere except in the interior

	updateAlpha();
	updateBeta();
	updateGamma();

	// Initialize the a coefficients. 9 Matrices to be created.
	//updateCoeff();

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


	// _a1 = _beta(i,j)/2/_deltaCsi/_deltaEta/-2.0/(_alpha(i,j)/_deltaCsi + _gamma(i,j)/_deltaEta);
	// _a2 = - _gamma(i,j)/pow2(_deltaEta)/-2/(_alpha(i,j)/_deltaCsi + _gamma(i,j)/_deltaEta);
	// _a4 = - _alpha(i,j)/pow2(_deltaEta)/-2/(_alpha(i,j)/_deltaCsi + _gamma(i,j)/_deltaEta);



}

void Problem::solve2()
{
	unsigned int iter = 0;
	while(iter < 20)
	{	
		for (int j = 1; j < _N + 1; ++j)
		{
			for (int i = 1; i < _M + 1; ++i)
			{
				updateAlpha();
				updateBeta();
				updateGamma();
				_a1 =  _beta(i,j)/2/_deltaCsi/_deltaEta/-2.0/(_alpha(i,j)/pow2(_deltaCsi) + _gamma(i,j)/pow2(_deltaEta));
				_a2 = -_gamma(i,j)/pow2(_deltaEta)/-2/(_alpha(i,j)/pow2(_deltaCsi) + _gamma(i,j)/pow2(_deltaEta));
				_a4 = -_alpha(i,j)/pow2(_deltaEta)/-2/(_alpha(i,j)/pow2(_deltaCsi) + _gamma(i,j)/pow2(_deltaEta));

				_x(i,j) = _a1 * ( _x(i+1,j+1) - _x(i-1,j+1) - _x(i+1,j-1) - _x(i-1, j-1) ) + _a2 * (_x(i,j+1) + _x(i,j-1)) + _a4 * (_x(i+1,j) + _x(i-1,j));
				std::cout << _x(i,j) << std::endl;
			}
		}

		for (int i = 1; i < _N + 1; ++i)
		{
			for (int j = 1; j < _M + 1; ++j)
			{
				updateAlpha();
				updateBeta();
				updateGamma();
				_a1 = _beta(i,j)/2/_deltaCsi/_deltaEta/-2.0/(_alpha(i,j)/pow2(_deltaCsi) + _gamma(i,j)/pow2(_deltaEta));
				_a2 = - _gamma(i,j)/pow2(_deltaEta)/-2/(_alpha(i,j)/pow2(_deltaCsi) + _gamma(i,j)/pow2(_deltaEta));
				_a4 = - _alpha(i,j)/pow2(_deltaEta)/-2/(_alpha(i,j)/pow2(_deltaCsi) + _gamma(i,j)/pow2(_deltaEta));

				_y(i,j) = _a1 * ( _y(i+1,j+1) - _y(i-1,j+1) - _y(i+1,j-1) - _y(i-1, j-1) ) + _a2 * (_y(i,j+1) + _y(i,j-1)) + _a4 * (_y(i+1,j) + _y(i-1,j));
				std::cout << _y(i,j) << std::endl;

			}
		}

		++iter;
	}	
	_x.print();
	_y.print();
}

void Problem::solve()
{
	unsigned int iterX= 0;

	while(iterX < 20)
	{
		sweepCols('x');
		++iterX;
		std::cout << "Iter = " << iterX << std::endl;
	}

	unsigned int iterY= 0;

	while(iterY < 20)
	{
		sweepCols('y');
		++iterY;
		std::cout << "Iter = " << iterY << std::endl;
	}


	// sweepCols();
}


void Problem::sweepCols(char var)
{

	// At each column of the domain (from j=1 to j=M_)
	// A matrix and a vector are created. Then the Gauss-Seidel algorithm is used.
	
	if (var == 'x')
	{
	
		std::cout << "Sweeping x over the cols " << std::endl;
		for (int j = _M ; j > 0; --j)   // sweep over all jnternal columns. From j=_M to j=1
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

			//updateCoeff();
		}
		updateAlpha();
		updateBeta();
		updateGamma();

		std::cout << "Sweeping x forward over the cols " << std::endl;
		for (int j = 1 ; j < _M + 1; ++j)   // sweep over all jnternal columns. From j=_M to j=1
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

			//updateCoeff();
		}
		updateAlpha();
		updateBeta();
		updateGamma();
	}

	if (var == 'y')
	{
		
		std::cout << "Sweeping y over the cols " << std::endl;
		for (int j = _M ; j > 0; --j)   // sweep over all jnternal rows. From j=1 to j=_M
		{
			Matrix<double> Ayj(_N, _N); // Renamed Ay into Aj to highlight the dependency on the column
			Vector<double> byj(_N);

			buildColMatrix(Ayj, j);
			buildColVector_y(byj, j);
			std::cout << "Ayj : " << std::endl;
			Ayj.print();

			std::cout << "byj : " << std::endl;
			byj.print();

			// Store the solution
			Vector<double> solyj = GaussSeidel(Ayj,byj);

			std::cout << "sol " << std::endl;
			solyj.print();

			// Insert solution inside the domain
			for (int i = 1; i < _N + 1; ++i)    // Insert solution
			{
				_y(i,j) = solyj(i-1);
			}


			
			std::cout << "y : " << std::endl;
			_y.print();

			//updateCoeff();
		}
		updateAlpha();
		updateBeta();
		updateGamma();

		std::cout << "Sweeping y forward over the cols " << std::endl;
		for (int j = 1 ; j < _M + 1; ++j)   // sweep over all jnternal rows. From j=1 to j=_M
		{
			Matrix<double> Ayj(_N, _N); // Renamed Ay into Aj to highlight the dependency on the column
			Vector<double> byj(_N);

			buildColMatrix(Ayj, j);
			buildColVector_y(byj, j);
			std::cout << "Ayj : " << std::endl;
			Ayj.print();

			std::cout << "byj : " << std::endl;
			byj.print();

			// Store the solution
			Vector<double> solyj = GaussSeidel(Ayj,byj);

			std::cout << "sol " << std::endl;
			solyj.print();

			// Insert solution inside the domain
			for (int i = 1; i < _N + 1; ++i)    // Insert solution
			{
				_y(i,j) = solyj(i-1);
			}


			
			std::cout << "y : " << std::endl;
			_y.print();

			//updateCoeff();
		}
		updateAlpha();
		updateBeta();
		updateGamma();
	}
}

void Problem::sweepRows(char var)
{
	// sweep over all internal rows. From i=1 to i=_N (i < _N + 1)
	// build a linear system in each row (matrix + vector)
	// solve the linear system and store the solution
	// Fill the domain with the solution
	if (var =='x')
	{

		std::cout << "Sweeping x over the rows " << std::endl;
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

			//updateCoeff();
		}
		updateAlpha();
		updateBeta();
		updateGamma();

		std::cout << "Sweeping x backwards over the rows " << std::endl;
		for (int i = _N; i > 0; --i)   
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

			//updateCoeff();
		}
		updateAlpha();
		updateBeta();
		updateGamma();
	}

	if (var =='y')
	{
	
		std::cout << "Sweeping y over the rows " << std::endl;

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

			//updateCoeff();
		}
		updateAlpha();
		updateBeta();
		updateGamma();

		for (int i = _N; i > 0; --i)
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

			//updateCoeff();
		}
		updateAlpha();
		updateBeta();
		updateGamma();

		

	}
}

void Problem::updateAlpha()
{
	std::cout << "Updating alpha " << std::endl;
	for (int i = 1; i < _alpha.rows() - 1; ++i)
	{
		for (int j = 1; j < _alpha.cols() - 1; ++j)
		{
			
			_alpha(i,j) = 1.0/4.0/pow2(_deltaEta) *( pow2( _x(i,j + 1) - _x(i,j - 1)) + pow2(_y(i,j + 1) - _y(i,j-1)) );
			
		}
	}
	std::cout << "alpha updated" << std::endl;
	_alpha.print();

}

void Problem::updateBeta()
{
std::cout << "Updating beta " << std::endl;
	for (int i = 1; i < _beta.rows() - 1; ++i)
	{
		for (int j = 1; j < _beta.cols() - 1; ++j)
		{	
			
			_beta(i,j) = 1.0/4./_deltaEta/_deltaCsi*( (_x(i+1,j) - _x(i-1,j))*(_x(i,j+1) - _x(i, j -1)) + (_y(i+1,j) - _y(i-1,j))*(_y(i,j+1) - _y(i,j-1)) );	
		}
	}
	std::cout << "Updated beta : " << std::endl;
	_beta.print();

}


void Problem::updateGamma()
{
	std::cout << "Updating gamma " << std::endl;
	for (int i = 1; i < _alpha.rows() - 1; ++i)
	{
		for (int j = 1; j < _alpha.cols() - 1; ++j)
		{
			_gamma(i,j) = 1.0/4.0/pow2(_deltaCsi) *(pow2( _x(i + 1,j ) - _x(i - 1,j )) + pow2(_y(i + 1,j) - _y(i - 1,j)));	
		}
	}

	std::cout << "Updated gamma : " << std::endl;
	 _gamma.print() ;
}


// void Problem::updateCoeff()
// {
// 	std::cout << "Updating coeff " << std::endl;
// 	_a1 = _beta(i,j)/2/_deltaCsi/_deltaEta/-2.0/(_alpha(i,j)/_deltaCsi + _gamma(i,j)/_deltaEta);
// 	_a2 = - _gamma(i,j)/pow2(_deltaEta)/-2/(_alpha(i,j)/_deltaCsi + _gamma(i,j)/_deltaEta);
// 	_a4 = - _alpha(i,j)/pow2(_deltaEta)/-2/(_alpha(i,j)/_deltaCsi + _gamma(i,j)/_deltaEta);
	
// }


void Problem::buildRowMatrix(Matrix<double>& A, unsigned int i)
{
	for (int l = 0; l < A.cols(); ++l)
	{
		int k = l + 1;
		A(l,l) = -2*_alpha(i,k)/pow2(_deltaCsi) - 2*_gamma(i,k)/pow2(_deltaCsi);   // set diagonal. There should be an index that travels internally to the domain, like k = l +1
		if (l < _N - 1)
		{
			A(l,l+1) = _gamma(i,k)/pow2(_deltaEta);
		}
		if (l > 0)
		{
			A(l,l-1) = _gamma(i,k)/pow2(_deltaEta);
		}
	}

}

void Problem::buildRowVector_x(Vector<double>& bx, unsigned int i)  // integer i locates the row of the domain
{
	for (unsigned int j = 1; j < bx.size() + 1; ++j) // loop over internal points of the domain. Given a j, l is determined as l = j - 1. l is used to fill the b vector
	{
		unsigned int l = j - 1;
		bx(l) =_beta(i,j)/2/_deltaCsi/_deltaEta*(_x(i+1, j+1) - _x(i+1, j-1) - _x(i-1, j+1) + _x(i-1, j-1) ) - _alpha(i,j)/pow2(_deltaCsi)*(_x(i+1,j) + _x(i-1,j)); 
	}
}

void Problem::buildRowVector_y(Vector<double>& by, unsigned int i)  // integer i locates the row of the domain
{
	for (unsigned int j = 1; j < by.size() + 1; ++j) // loop over internal points of the domain. Given a j, l is determined as l = j - 1. l is used to fill the b vector
	{
		unsigned int l = j - 1;
		std::cout << "_beta(i,j)/2/_deltaCsi/_deltaEta*(_y(i+1, j+1) - _y(i+1, j-1) - _y(i-1, j+1) + _y(i-1, j-1) )  =  " << _beta(i,j)/2/_deltaCsi/_deltaEta*(_y(i+1, j+1) - _y(i+1, j-1) - _y(i-1, j+1) + _y(i-1, j-1) ) << std::endl;
		std::cout << "_alpha(i,j)/pow2(_deltaEta)*(_y(i+1,j) + _y(i-1,j)) " << _alpha(i,j)/pow2(_deltaEta)*(_y(i+1,j) + _y(i-1,j)) << std::endl;
		by(l) = _beta(i,j)/2/_deltaCsi/_deltaEta*(_y(i+1, j+1) - _y(i+1, j-1) - _y(i-1, j+1) + _y(i-1, j-1) ) - _alpha(i,j)/pow2(_deltaEta)*(_y(i+1,j) + _y(i-1,j)); 
	}
}



void Problem::buildColMatrix(Matrix<double>& A, unsigned int j)
{
	for (int l = 0; l < A.rows(); ++l)
	{
		int k = l + 1;
		A(l,l) = -2*_alpha(k,j)/pow2(_deltaCsi) - 2*_gamma(k,j)/pow2(_deltaCsi);   // set diagonal. There should be an index that travels internally to the domain, like k = l +1
		if (l < _N - 1)
		{
			A(l,l+1) = _alpha(k,j)/pow2(_deltaCsi);  // modified
		}
		if (l > 0)
		{
			A(l,l-1) = _alpha(k,j)/pow2(_deltaCsi);
		}
	}

}

// temporary. I don't remember why.
void Problem::buildColVector_y(Vector<double>& by, unsigned int j)  // integer j locates the column of the domain
{
	for (unsigned int i = 1; i < by.size() + 1; ++i) // loop over internal points of the domain. Given a i, l is determined as l = i - 1. l is used to fill the b vector
	{
		unsigned int l = i - 1;
		// std::cout << "_aE(i,j) = " << _aE(i,j) << std::endl;
		// std::cout << "_x(i,j+1) = " << _x(i,j+1) << std::endl;
		// std::cout << "_x(i,j-1) = " << _x(i,j-1) << std::endl;
		// std::cout << "_aNE(i,j) = " << _aNE(i,j) << std::endl;
		// std::cout << "x(i+1,j+1) = " << _x(i+1,j+1) << std::endl;
		// std::cout << "x(i-1,j+1) = " << _x(i-1,j+1) << std::endl;
		// std::cout << "x(i+1,j-1) = " << _x(i+1,j-1) << std::endl;
		// std::cout << "x(i-1,j-1) = " << _x(i-1,j-1) << std::endl;
		// std::cout << "_aNE(i,j)*(_x(i+1, j+1) - _x(i+1, j-1) - _x(i-1, j+1) + _x(i-1, j-1) ) = " << _aNE(i,j)*(_x(i+1, j+1) - _x(i+1, j-1) - _x(i-1, j+1) + _x(i-1, j-1) ) << std::endl;
		// std::cout << "_aE(i,j)*(_x(i,j+1) + _x(i,j-1)) = " << _aE(i,j)*(_x(i,j+1) + _x(i,j-1)) << std::endl;
		by(l) = _beta(i,j)/2/_deltaCsi/_deltaEta*(_y(i+1, j+1) - _y(i-1, j+1) - _y(i+1, j-1) + _y(i-1, j-1) ) - _gamma(i,j)/pow2(_deltaEta)*(_y(i,j+1) + _y(i,j-1)); 
		
	}
}

void Problem::buildColVector_x(Vector<double>& bx, unsigned int j)  // integer j locates the column of the domain
{
	for (unsigned int i = 1; i < bx.size() + 1; ++i) // loop over internal points of the domain. Given a i, l is determined as l = i - 1. l is used to fill the b vector
	{
		unsigned int l = i - 1;
		// std::cout << "_aE(i,j) = " << _aE(i,j) << std::endl;
		// std::cout << "_x(i,j+1) = " << _x(i,j+1) << std::endl;
		// std::cout << "_x(i,j-1) = " << _x(i,j-1) << std::endl;
		// std::cout << "_aNE(i,j) = " << _aNE(i,j) << std::endl;
		// std::cout << "x(i+1,j+1) = " << _x(i+1,j+1) << std::endl;
		// std::cout << "x(i-1,j+1) = " << _x(i-1,j+1) << std::endl;
		// std::cout << "x(i+1,j-1) = " << _x(i+1,j-1) << std::endl;
		// std::cout << "x(i-1,j-1) = " << _x(i-1,j-1) << std::endl;
		// std::cout << "_aNE(i,j)*(_x(i+1, j+1) - _x(i+1, j-1) - _x(i-1, j+1) + _x(i-1, j-1) ) = " << _aNE(i,j)*(_x(i+1, j+1) - _x(i+1, j-1) - _x(i-1, j+1) + _x(i-1, j-1) ) << std::endl;
		// std::cout << "_aE(i,j)*(_x(i,j+1) + _x(i,j-1)) = " << _aE(i,j)*(_x(i,j+1) + _x(i,j-1)) << std::endl;
		bx(l) = _beta(i,j)/2/_deltaCsi/_deltaEta*(_x(i+1, j+1) - _x(i-1, j+1) - _x(i+1, j-1) + _x(i-1, j-1) ) - _gamma(i,j)/pow2(_deltaEta)*(_x(i,j+1) + _x(i,j-1)); 
		// std::cout << "bx(l) = " << bx(l) << std::endl;
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


void Problem::save()
{
	_y.save("y");
	_x.save("x");
}

