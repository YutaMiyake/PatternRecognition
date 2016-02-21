#ifndef __MATRIX__H__
#define __MATRIX__H__

#include <cassert>
#include <vector>

typedef std::vector< std::vector<double> > Matrix;

Matrix matrix_sum(const Matrix& A, const Matrix& B);
Matrix scalar_mul(const Matrix& A, double a);
Matrix matrix_inv(const Matrix& A);
double matrix_det(const Matrix& A);
std::vector<double> dot_mat_vec(Matrix A, std::vector<double> v);
double dot_vecs(std::vector<double> v,std::vector<double> u);

#endif

