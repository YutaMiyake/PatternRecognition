#ifndef __UTIL__H__
#define __UTIL__H__

#include <vector>
#include <cassert>

namespace util{
  std::vector<double> range(double a, double b, double c) {
      std::vector<double> array;
      while(a <= c) {
          array.push_back(a);
          a += b;
      }
      return array;
  }

  std::vector<double> linspace(double a, double b, int n)
  {
    assert(n >= 0); // non-negative

    std::vector<double> array;
    double step = (b-a) / (n-1);

    while(a <= b) {
        array.push_back(a);
        a += step;
    }
    return array;
  }
}

#endif