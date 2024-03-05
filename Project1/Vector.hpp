#ifndef VECTOR_H
#define VECTOR_H

#include <cassert>
#include <fstream>
#include <vector>
#include <iomanip>
#include <iostream>

using namespace std;

/**
 * 
 *    Vector class with custom access operators, printing and saving into files
 * 
 */
template <typename T>
class Vector
{
public:
  Vector(const unsigned int N) : v(N,0), N(N) {}
  Vector() {}

  const T & operator()(const unsigned int i) const
  {
    assert(i < N);
    return v[i];
  }
  T & operator()(const unsigned int i)
  {
    assert(i < N);
    return v[i];
  }
  const T & operator[](const unsigned int i) const
  {
    assert(i < N);
    return v[i];
  }
  T & operator[](const unsigned int i)
  {
    assert(i < N);
    return v[i];
  }

  unsigned int size() const
  {
    return N;
  }

  // Prints the vector
  void print(const string prefix = "", const bool newline = false, const unsigned int pr = 6) const
  {
    if (prefix.length() != 0)
      cout << prefix << endl;
    for (unsigned int i = 0; i < v.size(); ++i)
      cout << showpos << scientific << setprecision(pr) << v[i] << " ";
    cout << endl;
    if (newline)
      cout << endl;
  }

  Vector<T>& operator=(Vector<T>& v)
  {
    return v;
  }


private:
  vector<T> v;
  const unsigned int N = 0;
};

#endif /* VECTOR_H */