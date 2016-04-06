#include "Matrix.cpp"

int main(){
  double ** cov = new double*[4];

  for(int i = 0; i < 4; i++){
    cov[i] = new double[4];
  }
  cov[0][0] = 0;
  cov[0][1] = 0;
  cov[0][2] = 0;
  cov[0][3] = 0;

  cov[1][0] = 0;
  cov[1][1] = 1;
  cov[1][2] = 2;
  cov[1][3] = -1;

  cov[2][0] = 0;
  cov[2][1] = 2;
  cov[2][2] = -1;
  cov[2][3] = 4;

  cov[3][0] = 0;
  cov[3][1] = -1;
  cov[3][2] = 4;
  cov[3][3] = 2;

  print_matrix(subMatrix(toMatrix(cov,4,4),1,1,4,4));
  print_matrix(addZeroLoc(subMatrix(toMatrix(cov,4,4),1,1,4,4)));
}