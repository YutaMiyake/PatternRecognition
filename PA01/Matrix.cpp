#include "Matrix.h"

Matrix matrix_sum(const Matrix& A, const Matrix& B)
/* Pre: A and B are non-empty matrices with the same size. Return: A+B.*/
{
  assert(A.size() == B.size());
  assert(A[0].size() == B[0].size());

  int nrows = A.size();
  int ncols = A[0].size();
  Matrix matrix(nrows, std::vector<double>(ncols));
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
    matrix[i][j] = A[i][j] + B[i][j];
    }
  }
  return matrix;
}

Matrix scalar_mul(const Matrix& A, double a)
/*  Return: aA. */
{
  int nrows = A.size();
  int ncols = A[0].size();

  Matrix matrix(nrows, std::vector<double>(ncols));
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
    matrix[i][j] = a*A[i][j];
    }
  }
  return matrix;
}

Matrix matrix_inv(const Matrix& A){
  int nrows = A.size();
  int ncols = A[0].size();
  assert(nrows== 2);
  assert(ncols == 2);

  Matrix inv(nrows, std::vector<double>(ncols));
  double det = matrix_det(A);

  inv[0][0] = A[1][1]/det;
  inv[0][1] = -A[0][1]/det;
  inv[1][0] = -A[1][0]/det;
  inv[1][1] = A[0][0]/det;
  return inv;
}

double matrix_det(const Matrix& A){
  int nrows = A.size();
  int ncols = A[0].size();
  assert(nrows== 2);
  assert(ncols == 2);
  return A[0][0]*A[1][1] - A[0][1]*A[1][0];
}

std::vector<double> dot_mat_vec(Matrix A, std::vector<double> v){
  assert(A[0].size() == v.size());

  std::vector<double> dot(A.size());
  for(int row = 0; row < 2; row++){
    for(int col = 0; col < 2; col++){
      dot[row] += A[row][col]*v[col];
    }
  }
  return dot;
}
double dot_vecs(std::vector<double> v,std::vector<double> u){
  assert(v.size() == u.size());
  double dot = 0.0;

  for(int row = 0; row < v.size(); row++){
    dot += v[row]*u[row];
  }
  return dot;
}
