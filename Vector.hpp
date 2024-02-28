#ifndef VECTOR_H
#define VECTOR_H

#include <cassert>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

/**
 * Class that stores a 1D vector and enables printing and saving.
 */
template <typename T>
class Vector
{
public:
  Vector(const unsigned int N) : v(N), N(N) {}
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

  // // Saves the vector
  // void save(const string filename, const unsigned int pr = 12) const
  // {
  //   ofstream f;
  //   f.open(filename);
  //   for (unsigned int i = 0; i < v.size(); ++i)
  //     f << scientific << v[i] << endl;
  //   f.close();
  // }

private:
  vector<T> v;
  const unsigned int N = 0;
};

#endif /* VECTOR_H */