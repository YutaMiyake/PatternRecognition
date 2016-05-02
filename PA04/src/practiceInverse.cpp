#include "Matrix.cpp"
#include <iostream>

int main(){
  Matrix ex1 = {{2,4},{3,2}};
  Matrix ex2 = {{1,8,4},{1,7,54},{-1,2,-3}};
  Matrix ex3 =  {{7,4,2,10,6 },
                {-12,-4,6,3,8 },
                {1,0,4,2,1},
                {-28,4,2,71,-1},
                {910,71,4,7,1}};

  print_matrix(inverse(ex1));
  print_matrix(inverse(ex2));

  return 0;
}