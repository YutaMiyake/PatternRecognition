#include "Matrix.cpp"

int main(){
  Matrix m1(2,std::vector<double>(2));
  Matrix m2(2,std::vector<double>(2));

  m1[0][0] = 1;
  m1[0][1] = 2;
  m1[1][0] = 3;
  m1[1][1] = 4;

  m2[0][0] = 5;
  m2[0][1] = 6;
  m2[1][0] = 7;
  m2[1][1] = 8;

  print_matrix(matrix_mul_matrix(m1,m2));

  return 0;
}