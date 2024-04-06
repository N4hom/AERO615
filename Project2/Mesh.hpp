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

struct Face
{
	double delta_ = 0;
	
}; 

class Cell
{
	double xc_ = 0;
	double yc_ = 0;
	const unsigned int dim_;

	Vector<double> xFaces_;
	Vector<double> yFaces_;
public:
	Cell();
	~Cell(){};
	
};

Cell::Cell():
dim_(2),
xFaces_(dim_),
yFaces_(dim_)
{}

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


public:
	Mesh(unsigned int N, unsigned int M, std::string filename_x , std::string filename_y);
	~Mesh();

	friend std::vector<double> split(const std::string &s, char delimiter);
	bool read(std::string filename, Matrix<double>& xy);

	Matrix<double>& x();
	Matrix<double>& y();
	Matrix<double>& xC();
	Matrix<double>& yC();
};

Mesh::Mesh(unsigned int N, unsigned int M, std::string filename_x, std::string filename_y):
x_(N,M),
y_(N,M),
xC_(N+1,M+1),
yC_(N+1,M+1),
Imax_(N-1),
Jmax_(M-1),
Icmax_(N),
Jcmax_(M),
N_(N),
M_(M),
Nc_(N + 1),
Mc_(M + 1)
{
	read(filename_x, x_);
	read(filename_y, y_);

	for (unsigned int i = 0; i < x_.rows(); ++i)
	{
		for (unsigned int j = 0; j < x_.cols(); ++j)
		{
			if (j > 0)
			{
				/* code */
				xC_(i,j) = 0.5*(x_(i,j) + x_(i,j-1));
			}

			// lower boundary
			if (i == x_.rows() - 1)
			{
				xC_(i + 1 , j) = xC_(i,j);
			}

			// right boundary
			if (j == x_.cols() - 1)
			{
				xC_(i,j+1) = x_(i,j);

				// bottom-right corner
				if (i == x_.rows() - 1)
				{
					xC_(i+1,j+1) = x_(i,j);
				}
			}

		}
	}

	xC_.print();



	for (unsigned int i = 0; i < N + 1; ++i)
	{
		for (int j = 0; j < M + 1; ++j)
		{
			if (i > 0 && j > 0 && i < Icmax_ && j < Jcmax_)
			{
				yC_(i,j) = 0.25*(y_(i-1,j-1) + y_(i,j-1) + y_(i-1,j) + y_(i,j));
				if (j == 0)
				{
					yC_(i,j) = 0.5*(y_(i+1,j) + y_(i,j));
				}			
			}
			else if (i == 0 )
			{
				if (j == 0)
				{
					yC_(i,j) = y_(i,j);
				}
				else if (j == Jcmax_)
				{
					yC_(i,j) = y_(i,Jmax_);
				}
				else
				{
					yC_(i,j) = 0.5*(y_(i,j) + y_(i,j-1));
				}
			}
			else if (i == Icmax_)
			{
				if (j == 0)
				{
					yC_(i,j) = y_(Imax_,j);
				}
				else if (j == Jcmax_)
				{
					yC_(i,j) = y_(Imax_,Jmax_);
				}
				else
				{
					yC_(i,j) = 0.5*(y_(Imax_,j) + y_(Imax_,j-1));
				}
			}
			else if (j == 0)
			{
				if (i == 0)
				{
					// Already set
				}
				else if (i == Icmax_)
				{
					yC_(i,j) = y_(Imax_,j);
				}
				else
				{
					yC_(i,j) = 0.5*(y_(i,0) + y_(i-1,0));
				}
			}
			else if (j == Jcmax_)
			{
				if (i == Icmax_)
				{
					// Already set
				}
				else if (i == Icmax_)
				{
					//
				}
				else
				{
					yC_(i,j) = 0.5*(y_(i,0) + y_(i-1,0));
				}
			}
		}
	}

	yC_.print();
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
