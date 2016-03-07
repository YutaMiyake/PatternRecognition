#include "rgb.h"
RGB::RGB(){}
RGB::RGB(int r, int g, int b){
  this->r = r;
  this->g = g;
  this->b = b;
}
RGB& RGB::operator=(const RGB &val){
  if(this == &val){
    return *this;
  }
  this->r = val.r;
  this->g = val.g;
  this->b = val.b;
  return *this;
}
bool RGB::operator==(const RGB &val){
  return this->r == val.r && this->g == val.g && this->b == val.b;
}
bool RGB::operator!=(const RGB &val){
  return this->r != val.r || this->g != val.g || this->b != val.b;
}
std::ostream& operator<<(std::ostream& os, const RGB &val){
  os << "(" << val.r << "," << val.g << "," << val.b << ")";
  return os;
}