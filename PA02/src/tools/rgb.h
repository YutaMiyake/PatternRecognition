#ifndef RGB_H
#define RGB_H

#include <iostream>

struct RGB {
  RGB();
  RGB(int, int, int);
  RGB& operator=(const RGB &val);
  bool operator==(const RGB &val);
  bool operator!=(const RGB &val);
  friend std::ostream& operator<<(std::ostream &os, const RGB &val);
  int r, g, b;
};

#endif
