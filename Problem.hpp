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

	//- Lenght and height
	double _L = 5.0, _H = 1.0;   
	
	//- Number of unknowns along each direction per variable
	unsigned int _N, _M;  
	unsigned int _Imax; 
	unsigned int _Jmax;
	
	//- Grid size along csi and eta
	double _deltaCsi, _deltaEta;
	
	//- Solution storage
	Matrix<double> _x, _y;  
	
	// Coefficient storage at location ij
	double _alpha, _beta, _gamma;  
	double _a1, _a2, _a4;

	// Tolerance for convergence criterion
	double _tol = 1e-5;   

public:
	Problem(unsigned int, unsigned int);
	~Problem();

	void initialize();
	void updateCoeff();
	void updateAlpha(unsigned int i, unsigned int j);
	void updateBeta(unsigned int i, unsigned int j);
	void updateGamma(unsigned int i, unsigned int j);
	void solve();
	void save();

	Matrix<double>& x(){ return _x;};
	Matrix<double>& y(){ return _y;};

	friend double pow2(double f);
};

#endif 