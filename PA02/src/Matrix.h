#ifndef __MATRIX__H__
#define __MATRIX__H__

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

typedef std::vector< std::vector<double> > Matrix;

Matrix matrix_sum(const Matrix& A, const Matrix& B);
Matrix scalar_mul(const Matrix& A, double a);
Matrix matrix_inv(const Matrix& A);
double matrix_det(const Matrix& A);
std::vector<double> dot_mat_vec(Matrix A, std::vector<double> v);
std::vector<double> vec_sum(std::vector<double> v1, std::vector<double> v2);
std::vector<double> vec_diff(std::vector<double> v1, std::vector<double> v2);
double dot_vecs(std::vector<double> v,std::vector<double> u);
void loadFile(Matrix &points, std::string filename);
void print_vec(std::vector<double> vec);
void print_matrix(Matrix mat);
Matrix randSample(Matrix data, int num);

#endif

