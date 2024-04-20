#ifndef MESH_
#define MESH_

#include "Matrix.hpp"
#include "Vector.hpp"

std::vector<double> split(const std::string &s, char delimiter) {
    std::vector<double> tokens;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(std::stod(token));
    }
    return tokens;
}

// Forward declaration
class Problem;
class Variable;

class Mesh
{
private:

	Matrix<double> x_;
	Matrix<double> y_;

	Matrix<double> xC_;
	Matrix<double> yC_;

	unsigned int Jmax_;
	unsigned int Imax_;
	unsigned int Jcmax_;
	unsigned int Icmax_;

	unsigned int N_;
	unsigned int M_;

	unsigned int Nc_;
	unsigned int Mc_;

	// Matrix<double> xFacesLeft_;
	// Matrix<double> xFacesRight_;
	// Matrix<double> xFacesTop_;
	// Matrix<double> xFacesBottom_;

	Matrix<Coefficients> xFaces_;
	Matrix<Coefficients> yFaces_;

	// Matrix<double> yFacesLeft_;
	// Matrix<double> yFacesRight_;
	// Matrix<double> yFacesTop_;
	// Matrix<double> yFacesBottom_;

	Matrix<double> area_;


	friend class Problem;
	friend class Variable;
public:
	Mesh(unsigned int N, unsigned int M, std::string filename_x , std::string filename_y);
	~Mesh();

	friend std::vector<double> split(const std::string &s, char delimiter);
	bool read(std::string filename, Matrix<double>& xy);

	Matrix<double>& x();
	Matrix<double>& y();
	Matrix<double>& xC();
	Matrix<double>& yC();
	Matrix<double> area();

};

//
//             ____ ____ ____ ____
//            |	   |	|	 |	  |
//            |____|____|____|____|
//            |    |    |    |    |
//            |____|____|____|____|
//            |    |    |    |    |
//            |____|____|____|____|
//
//
//
Mesh::Mesh(unsigned int N, unsigned int M, std::string filename_x, std::string filename_y):
x_(N,M),
y_(N,M),
xC_(N+3,M+3),   // NUMBER OF CELLS = NUMBER OF NODES - 1 + NUMBER OF GHOST CELLS(4)
yC_(N+3,M+3),
Imax_(N-1),
Jmax_(M-1),
Icmax_(N+2),    // CELL INDEX MAX = NUMBER OF CELLS - 1
Jcmax_(M+2),
N_(N),
M_(M),
Nc_(N - 1),
Mc_(M - 1),
xFaces_(Nc_,Mc_),
yFaces_(Nc_,Mc_),
area_(N+3,M+3)
{

	std::cout << " Initializing geometry " << std::endl;

	read(filename_x, x_);
	read(filename_y, y_);

	for (unsigned int i = 0; i < Imax_ ; ++i)
	{
		for (unsigned int j = 0; j < Jmax_ ; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;

			xC_(ic,jc) = 0.25*(x_(i,j) + x_(i+1,j) + x_(i+1,j+1) + x_(i,j+1));

			if (j == 0)
			{
				xC_(ic,1) = - xC_(ic,2);
				xC_(ic,0) = - 2*xC_(ic,2);
			}

			if (i == 0)
			{
				xC_(1,jc) = xC_(2,jc);
				xC_(0,jc) = xC_(2,jc);
			}

			if (i == Imax_ - 1)
			{
				xC_(Icmax_ - 1,jc) = xC_(Icmax_ - 2, jc); 
				xC_(Icmax_ ,jc) = xC_(Icmax_ - 1, jc); 
			}

			if (j == Jmax_ - 1)
			{
				double deltaX = xC_(ic, jc ) - xC_(ic,jc-1); 

				xC_(ic,Jcmax_ - 1) = xC_(ic, Jcmax_-2) + deltaX;
				xC_(ic,Jcmax_    ) = xC_(ic, Jcmax_-1) + deltaX; 
			}

			
				// double deltaX = xC_(ic , Jcmax_ - 3) - xC_(ic , Jcmax_ - 2);
				// xC_(ic, Jcmax_ - 1) = xC_(ic, Jcmax_ - 2);
				// xC_(ic, Jcmax_) = xC_(ic, Jcmax_ - 1) ;

			// std::cout << "x(i,j) " << x_(i,j) << endl;
			// std::cout << "xc(i,j) " << xC_(ic,jc) << endl;
		}		
	}

	xC_.print();



	for (unsigned int i = 0; i < Imax_; ++i)
	{
		for (int j = 0; j < Jmax_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;

			yC_(ic,jc) = 0.25*(y_(i,j) + y_(i+1,j) + y_(i+1,j+1) + y_(i,j+1));

			if (j == 0)
			{
				yC_(ic,1) = yC_(ic,2);
				yC_(ic,0) = yC_(ic,1);
			}

			if (i == 0)
			{
				double deltaY = yC_(3,jc) - yC_(2,jc);
				

				yC_(1 , jc) = yC_(2,jc) - deltaY;
				yC_(0 , jc) = yC_(1,jc) - deltaY;
			}

			if (i == Imax_ - 1)
			{
				double deltaY = yC_(Icmax_ - 3 , jc ) - yC_(Icmax_ - 2 , jc );
				yC_(Icmax_ - 1 , jc) = yC_(Icmax_ - 2 , jc) - deltaY; 
				yC_(Icmax_,jc) = yC_(Icmax_ - 1 , jc) - deltaY; 
			}

			if (j == Jmax_ - 1)
			{

				yC_(ic,Jcmax_ - 1) = yC_(ic, Jcmax_-2) ;
				yC_(ic,Jcmax_    ) = yC_(ic, Jcmax_-1) ; 
			}


			
		}
	}

	for (unsigned int i = 0; i < Imax_; ++i)
	{
		for (unsigned int j = 0; j < Jmax_; ++j)
		{
			unsigned int ic = i + 2;
			unsigned int jc = j + 2;
			if (i == 0)
			{
				double deltaY = yC_(2,jc) - yC_(3,jc);
	
				yC_(1 , jc) = yC_(2,jc) + deltaY;
				yC_(0 , jc) = yC_(1,jc) + deltaY;


			}
		}

		yC_.print();
	}

	// Define faces

	// Store face areas
	for (unsigned int i = 0; i < Nc_ ; ++i)
	{
		for (unsigned int j = 0; j < Mc_ ; ++j)
		{
			xFaces_(i,j).w  =  y_(i + 1 , j     ) - y_(i     , j     );
			xFaces_(i,j).e  =  y_(i     , j + 1 ) - y_(i + 1 , j + 1 );
			xFaces_(i,j).n  =  y_(i     , j     ) - y_(i     , j + 1 );
			xFaces_(i,j).s  =  y_(i + 1 , j + 1 ) - y_(i + 1 , j     );

		}	
	}
	
	// Define left and right faces length component in the y direction
	for (unsigned int i = 0; i < Imax_ ; ++i)
	{
		for (unsigned int j = 0; j < Jmax_ ; ++j)
		{
			
			yFaces_(i,j).w  =  x_(i + 1 , j     ) - x_(i     , j     );
			yFaces_(i,j).e  =  x_(i     , j + 1 ) - x_(i + 1 , j + 1 );
			yFaces_(i,j).n  =  x_(i     , j     ) - x_(i     , j + 1 );
			yFaces_(i,j).s  =  x_(i + 1 , j + 1 ) - x_(i + 1 , j     );

		}	
	}		


	for (unsigned int i = 0; i < Imax_ ; ++i)
	{
		for (unsigned int j = 0; j < Jmax_ ; ++j)
		{
			std::cout << "i " << i << std::endl;
			std::cout << "j " << j << std::endl;
			std::cout << "xFaces_ " << "\n" << std::endl;
			xFaces_(i,j).print();
			std::cout << "yFaces_ " << "\n" << std::endl;
			yFaces_(i,j).print();
		}	
	}		

	
	// std::cout << "xFacesLeft_ " << "\n" << std::endl;
	// xFacesLeft_.print();
	// std::cout << "xFacesRight_ " << "\n" << std::endl;
	// xFacesRight_.print();
	// std::cout << "xFacesTop_ " << "\n" << std::endl;
	// xFacesTop_.print();
	// std::cout << "xFacesBottom_ " << "\n" << std::endl;
	// xFacesBottom_.print();

	// std::cout << "yFacesLeft_ " << "\n" << std::endl;
	// yFacesLeft_.print();
	// std::cout << "yFacesRight_ " << "\n" << std::endl;
	// yFacesRight_.print();
	// std::cout << "yFacesTop_ " << "\n" << std::endl;
	// yFacesTop_.print();
	// std::cout << "yFacesBottom_ " << "\n" << std::endl;
	// yFacesBottom_.print();

	double areaTot = 0;

	for (int i = 0; i < Nc_ ; ++i)
	{
		for (int j = 0; j < Mc_ ; ++j)
		{
			int ic = i + 2;
			int jc = j + 2;
			double Aij;

			Aij = 0.5 * ( (x_(i + 1 , j + 1) - x_(i  , j )) * (y_(i  , j + 1) - y_(i + 1, j )) 
				        - (y_(i + 1 , j + 1) - y_(i  , j )) * (x_(i  , j + 1) - x_(i + 1, j )) );
			

			area_(ic , jc) = Aij;

			areaTot += Aij;
		}
	}

	
	area_.print();

	std::cout << "Geometry initialized " << std::endl;
}

Mesh::~Mesh(){}

//- Get filename and matrix storing either x and y
//  and read the file to fill the matrix
bool Mesh::read(std::string filename, Matrix<double>& xy)
{
	std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file\n";
        return false;
    }

    std::string line;
    int i = 0;
    while (std::getline(file, line)) {
        // Split the line into values
        std::vector<double> row = split(line, ',');
        // Add the row to the matrix
        xy.setRow(i,row);
        ++i;
    }

    // Close the file
    file.close();

    return true;
}



Matrix<double>& Mesh::x()
{
	return x_;
}

Matrix<double>& Mesh::y()
{
	return y_;
}

Matrix<double>& Mesh::xC()
{
	return xC_;
}

Matrix<double>& Mesh::yC()
{
	return yC_;
}

#endif