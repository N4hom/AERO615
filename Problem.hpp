#ifndef PROBLEM
#define PROBLEM

#include <iostream>
#include <math.h>
#include <cmath>
#include "Matrix.hpp"
#include "Vector.hpp"

class Problem
{
private:

	double _L = 5.0, _H = 1.0;   // Lenght and height
	unsigned int _N, _M;  // Size of the problem
	unsigned int _Imax;
	unsigned int _Jmax;
	double _deltaCsi, _deltaEta;
	Matrix<double> _x, _y;  // Solution storage
	double _alpha, _beta, _gamma;  // Coefficient storage. They could be a vector but for now I'll keep them Matrices. I was kidding, they're matrices
	double _aP, _aEW, _aNS,_aC;   // Storage for linear system coefficient

	double _a1, _a2, _a4;
	double _tol = 1e-1;

	Matrix<double> _A; // Problem matrix


public:
	Problem(unsigned int, unsigned int);
	~Problem();
	void initialize();
	void updateCoeff();
	void updateAlpha(unsigned int i, unsigned int j);
	void updateBeta(unsigned int i, unsigned int j);
	void updateGamma(unsigned int i, unsigned int j);
	
	// Create linear system in each row of the domain. To be used inside sweepRows() where a temporary matrix is created

	Vector<double> GaussSeidel(Matrix<double>& A, Vector<double>& b, double tol=1e-7 , unsigned int N=5000);
	void solve2();

	// Matrix<double>& alpha(){ return _alpha;};
	// Matrix<double>& beta(){ return _beta;};
	// Matrix<double>& gamma(){ return _gamma;};

	Matrix<double>& x(){ return _x;};
	Matrix<double>& y(){ return _y;};

	friend double pow2(double f);

	void save();

	
};

#endif 