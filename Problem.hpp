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

	double _L = 1.0, _H = 1.0;   // Lenght and height
	unsigned int _N, _M;  // Size of the problem
	unsigned int _Imax;
	unsigned int _Jmax;
	double _deltaCsi, _deltaEta;
	Matrix<double> _x, _y;  // Solution storage
	Matrix<double> _alpha, _beta, _gamma;  // Coefficient storage. They could be a vector but for now I'll keep them Matrices. I was kidding, they're matrices
	double _aP, _aEW, _aNS,_aC;   // Storage for linear system coefficient

	double _a1, _a2, _a4;

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
	void sweepCols(char var);

	Vector<double> GaussSeidel(Matrix<double>& A, Vector<double>& b, double tol=1e-7 , unsigned int N=5000);
	void solve();
	void solve2();

	Matrix<double>& alpha(){ return _alpha;};
	Matrix<double>& beta(){ return _beta;};
	Matrix<double>& gamma(){ return _gamma;};

	Matrix<double>& x(){ return _x;};
	Matrix<double>& y(){ return _y;};

	friend double pow2(double f);

	void save();

	
};

#endif 