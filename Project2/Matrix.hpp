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
  Matrix(const unsigned int N, const unsigned int M) : N(N), M(M), A(N, vector<T>(M)) {}
        


  const unsigned int cols()
  {
     return M;
  }

  const unsigned int rows()
  {
     return N;
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

  // Assignment operator to copy values from another matrix
  Matrix<T>& operator=(const Matrix<T>& other) {
    if (this != &other) { // Check for self-assignment
        // Resize if necessary
        // Copy values
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                A[i][j] = other(i, j);
            }
        }
    }
    return *this;
}

// Addition operator to add two matrices
Matrix<T> operator+(const Matrix<T>& other) const {
    // Check if matrices have the same dimensions
    assert(N == other.rows() && M == other.cols());
    
    // Create a new matrix to store the result
    Matrix<T> result(N, M);
    
    // Perform element-wise addition
    for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < M; ++j) {
            result(i, j) = (*this)(i, j) + other(i, j);
        }
    }
    
    return result;
}

  // Print matrix
  void print() const
  {
    for (unsigned int i = 0; i < N; ++i)
    {
      for (unsigned int j = 0; j < M; ++j)
      {
          std::cout << A[i][j] << " " ;
      }

      std::cout << std::endl;
    }

    std::cout <<"\n" << std::endl;
    //std::cout << std::endl;
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

#endif 