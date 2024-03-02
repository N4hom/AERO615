#ifndef MATRIX_CLASS
#define MATRIX_CLASS

#define NDEBUG
#include <cassert>
#include <vector>
#include "Vector.hpp"
#include <iostream>

/**
 * Class that holds a N x M matrix with common matrix operations.
 */
template <typename T>
class Matrix {

private:
  // The size of this matrix
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

  // unsigned int rowNumber()
  // {
  //    return A.size();
  // }

  // unsigned int colNumber()
  // {
  //    return A[0].size();
  // }

  // Const access to the ith,jth element
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
    for (unsigned int j = 0; j < M; ++j)
    {
      for (unsigned int i = 0; i < N; ++i)
      {
        if (i > 0)
          f << ",";
        f << setprecision(pr) << A[i][j];
      }
      f << endl;
    }
    f.close();
  }


  // Gauss-Seidel algorithm

};

#endif 