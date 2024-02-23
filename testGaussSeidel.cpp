#include <iostream>
#include "Vector.hpp"
#include "Matrix.hpp"

// x = [0, 0, 0]                        
// a = [[4, 1, 2],[3, 5, 1],[1, 1, 3]]
// b = [4,7,3]
// solution = [0.5, 1.0, 0.5]

Vector<double> GaussSeidel(Matrix<double>& A, Vector<double>& b, double tol=1e-5 , unsigned int N=10)
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
		}

		++n;

	}

	return u;
}
int main(int argc, char const *argv[])
{
	Vector<double> b(3);
	b(0) = 4;
	b(1) = 7;
	b(2) = 3;

	std::vector<double> a0 = {4,1,2};
	std::vector<double> a1 = {3,5,1};
	std::vector<double> a2 = {1,1,3};

	Matrix<double> A(3,3);
	A.setRow(0,a0);
	A.setRow(1,a1);
	A.setRow(2,a2);

	A.print();

	Vector<double> u = GaussSeidel(A, b);
	u.print();

	return 0;
}