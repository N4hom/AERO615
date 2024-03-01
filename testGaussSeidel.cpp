#include <iostream>
#include "Vector.hpp"
#include "Matrix.hpp"

// x = [0, 0, 0]                        
// a = [[4, 1, 2],[3, 5, 1],[1, 1, 3]]
// b = [4,7,3]
// solution = [0.5, 1.0, 0.5]



// [10.0, -1.0, 2.0, 0.0],
//     [-1.0, 11.0, -1.0, 3.0],
//     [2.0, -1.0, 10.0, -1.0],
//     [0.0, 3.0, -1.0, 8.0],


// b = [6.0, 25.0, -11.0, 15.0]
Vector<double> GaussSeidel(Matrix<double>& A, Vector<double>& b, double tol = 1e-5, unsigned int N = 100)
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
    }

    return u;
}

template <typename T>
Vector<T> solveTridiagonal(Matrix<T>& A, Vector<T>& d) {
    unsigned int N = A.rows();
    Vector<T> x(N);

    // Ensure the matrix and vector sizes are compatible
    assert(N == A.cols() && N == d.size());

    // Create temporary vectors for the algorithm
    Vector<T> cPrime(N - 1);
    Vector<T> dPrime(N);

    // Forward sweep
    cPrime(0) = A(0, 1) / A(0, 0);
    dPrime(0) = d(0) / A(0, 0);

    for (unsigned int i = 1; i < N - 1; ++i) {
        T temp = A(i, i) - A(i, i - 1) * cPrime(i - 1);
        cPrime(i) = A(i, i + 1) / temp;
        dPrime(i) = (d(i) - A(i, i - 1) * dPrime(i - 1)) / temp;
    }

    // Backward substitution
    x(N - 1) = (d(N - 1) - A(N - 1, N - 2) * dPrime(N - 2)) / (A(N - 1, N - 1) - A(N - 1, N - 2) * cPrime(N - 2));

    for (int i = N - 2; i >= 0; --i) {
        x(i) = dPrime(i) - cPrime(i) * x(i + 1);
    }

    return x;
}

int main(int argc, char const *argv[])
{


	std::cout << "First test " << std::endl;
	
	std::vector<double> a0 = {4,1,2};
	std::vector<double> a1 = {3,5,1};
	std::vector<double> a2 = {1,1,3};

	Vector<double> b1(3);
	b1(0) = 4;
	b1(1) = 7;
	b1(2) = 3;


	Matrix<double> A1(3,3);
	A1.setRow(0,a0);
	A1.setRow(1,a1);
	A1.setRow(2,a2);

	A1.print();

	Vector<double> uTDMA = solveTridiagonal<double>(A1,b1);

	uTDMA.print();

	Vector<double> u1 = GaussSeidel(A1, b1);
	u1.print();


	std::cout << "second test " << std::endl;
	// [ 1.  2. -1.  1.]

	Vector<double> b2(4); 
	b2(0) = 6.0;
	b2(1) = 25.0;
	b2(2) = -11.0;
	b2(3) = 15.0;

	std::vector<double> a0_2 = {10.0, -1.0, 2.0, 0.0};
	std::vector<double> a1_2 = {-1.0, 11.0, -1.0, 3.0};
	std::vector<double> a2_2 = {2.0, -1.0, 10.0, -1.0};
	std::vector<double> a3_2 = {0.0, 3.0, -1.0, 8.0};

	Matrix<double> A2(4,4);
	A2.setRow(0,a0_2);
	A2.setRow(1,a1_2);
	A2.setRow(2,a2_2);
	A2.setRow(3,a3_2);

	A2.print();

	Vector<double> u2 = GaussSeidel(A2, b2);
	u2.print();


	std::vector<double> ay0 = {-6.125000e+02, +0.000000e+00, +0.000000e+00, +0.000000e+00, +0.000000e+00};
	std::vector<double> ay1 = {+0.000000e+00, -0.000000e+00, +0.000000e+00, +0.000000e+00, +0.000000e+00};
	std::vector<double> ay2 = {+0.000000e+00, -0.000000e+00, +0.000000e+00, +0.000000e+00, +0.000000e+00};
	std::vector<double> ay3 = {+0.000000e+00, -0.000000e+00, +0.000000e+00, +0.000000e+00, +0.000000e+00};
	std::vector<double> ay4 = {+0.000000e+00, +0.000000e+00, +0.000000e+00, +0.000000e+00, -1.813000e+03};

	


	Matrix<double> Ay(5,5);
	Ay.setRow(0,ay0);
	Ay.setRow(1,ay1);
	Ay.setRow(2,ay2);
	Ay.setRow(3,ay3);
	Ay.setRow(4,ay4);

	Vector<double> by(5);
	by(0) = -3.350735e+02;
	by(1) = -0.0;
	by(2) = -0.0;
	by(3) = -0.0;
	by(4) = -1.305476e+03;

    Vector<double> solution = GaussSeidel(Ay, by);

    // Print the solution
    solution.print() ;


	return 0;
}