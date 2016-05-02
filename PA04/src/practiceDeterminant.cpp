#include "Matrix.cpp"
#include <iostream>

int main(){
  Matrix ex1 = {{4,7},{8,9}};
  Matrix ex2 = {{7,6,8},{4,8,4},{7,5,3}};
  Matrix ex3 =  {{7,4,2,10,6 },
                {-12,-4,6,3,8 },
                {1,0,4,2,1},
                {-28,4,2,71,-1},
                {910,71,4,7,1}};

  std::cout << determinant(ex1) << std::endl;
  std::cout << determinant(ex2) << std::endl;
  std::cout << determinant(ex3) << std::endl;

  return 0;
}