#ifndef MATRIX_CLASS
#define MATRIX_CLASS

#define NDEBUG
#include <cassert>
#include <vector>
#include "Vector.hpp"
#include <iostream>

/**
 *
 *  Matrix class for data storage.
 * 
 */
template <typename T>
class Matrix {

private:
  //- size of this matrix
  const unsigned int N, M;

  // Matrix storage
  std::vector<std::vector<T> > A;


public:
  Matrix(unsigned int N, unsigned int M, T v = 0): 
        N(N), 
        M(M), 
        A(N, std::vector<T>(M, v))
        {}


  unsigned int cols()
  {
     return M;
  }

  unsigned int rows()
  {
     return N;
  }

  const unsigned int rows() const {
    return N;
  }

  const unsigned int cols() const {
    return M;
  }



  // Const access to the ith,jth element  (read only)
  const T &operator()(unsigned int i, unsigned int j) const
  {
    // Check if the indexes are consistent with the size
    assert(i < N && j < M);
    return A[i][j];
  }

  // Non const access for getting the ith, jth element
  T &operator()(unsigned int i, unsigned int j) 
  {

    // Check if the indexes are consistent with the size
    assert(i < N && j < M);
    return A[i][j];
  }

  // Operator for setting the entire matrix to a value
  void operator=(T v) 
  {
    for (unsigned int j = 0; j < M; ++j)
      setRow(j, v);
  }

  // Set the j-th row to v
  void setColumn(unsigned int j, T v)
  {
    assert(j < M);
    for (unsigned int i = 0; i < N; ++i)
      A[i][j] = v;
  }
  // Set the i-th column to v
  void setRow(unsigned int i, T v) 
  {
    assert(i < N);
    for (unsigned int j = 0; j < M; ++j)
      A[i][j] = v;
  }

  // Set the j-th row to vector v
  void setColumn(unsigned int j, std::vector<T> &v)
  {
    assert(j < M && vs.size() == N);
    for (unsigned int i = 0; i < N; ++i)
      A[i][j] = v[i];
  }

  // Set the j-th row to vector v
  void setColumn(unsigned int j, Vector<T> v)
  {
    assert(j < M && vs.size() == N);
    for (unsigned int i = 0; i < N; ++i)
      A[i][j] = v(i);
  }

  // Set the i-th column to vector v
  void setRow(unsigned int i, std::vector<T> &v) 
  {
    assert(i < N && vs.size() == M);
    for (unsigned int j = 0; j < M; ++j)
      A[i][j] = v[j];
  }

   // Set the i-th column to vector v
  void setRow(unsigned int i, Vector<T> v) 
  {
    assert(i < N && vs.size() == M);
    for (unsigned int j = 0; j < M; ++j)
      A[i][j] = v(j);
  }

  // Print matrix
  void print()
  {
    for (unsigned int i = 0; i < N; ++i)
    {
      for (unsigned int j = 0; j < M; ++j)
      {
          std::cout << A[i][j] << " " ;
      }

      std::cout << std::endl;
    }

    std::cout << std::endl;
  }

  // Saves the matrix in csv format
  void save(const string filename, const unsigned int pr = 12) const
  {
    ofstream f;
    f.open(filename);
    for (unsigned int i = 0; i < N; ++i)
    {
      for (unsigned int j = 0; j < M; ++j)
      {
        if (j > 0)
          f << ",";
        f << setprecision(pr) << A[i][j];
      }
      f << endl;
    }
    f.close();
  }



};


template <typename T>
void writeVTK(const Matrix<T>& x_matrix, const Matrix<T>& y_matrix, const string& filename) {
    assert(x_matrix.rows() == y_matrix.rows() && x_matrix.cols() == y_matrix.cols());

    unsigned int rows = x_matrix.rows();
    unsigned int cols = y_matrix.cols();

    // Open the file
    ofstream vtkfile(filename);
    if (!vtkfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    // Write the header
    vtkfile << "# vtk DataFile Version 3.0" << endl;
    vtkfile << "VTK output" << endl;
    vtkfile << "ASCII" << endl;
    vtkfile << "DATASET STRUCTURED_GRID" << endl;
    vtkfile << "DIMENSIONS   " << cols << "  " << rows << " 1 " << endl;
    vtkfile << "POINTS " << rows * cols << " FLOAT" << endl;

    

     // Write the points from bottom row to top row (reverse the row order)
    for (int i = rows - 1; i >= 0; --i) { // Start from the last row and go to the first row
        for (int j = 0; j < cols; ++j) { // Iterate columns left to right
            vtkfile << std::fixed << std::setprecision(14)
                    << x_matrix(i, j) << " " << y_matrix(i, j) << " 0.0" << std::endl;
        }
    }

    
    // Close the file
    vtkfile.close();
}

#endif 