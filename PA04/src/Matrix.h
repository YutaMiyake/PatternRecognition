#ifndef __MATRIX__H__
#define __MATRIX__H__

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

typedef std::vector< std::vector<double> > Matrix;

Matrix operator+(const Matrix& A, const Matrix& B);
Matrix operator*(const Matrix& A, double a);
Matrix operator/(const Matrix& A, double a);
Matrix matrix_mul_matrix(const Matrix& A,const Matrix& B);
Matrix inverse(const Matrix& A);
Matrix matrix_vec_diff(const Matrix& A, std::vector<double> v);
Matrix toMatrix(double** A, int rows, int cols);
Matrix transpose(const Matrix& A);
double** toPtr(const Matrix& A);

double determinant(const Matrix& A);
std::vector<double> dot_mat_vec(Matrix A, std::vector<double> v);
std::vector<double> operator+(std::vector<double> v1, std::vector<double> v2);
std::vector<double> operator-(std::vector<double> v1, std::vector<double> v2);
std::vector<double> operator*(std::vector<double> v1, double a);

std::vector<double> normalize(std::vector<double> v1);
std::vector<double> toVector(double* A, int size);
double dot_vecs(std::vector<double> v, std::vector<double> u);
double euclideanNorm(std::vector<double> v1);
double euclideanDistance(std::vector<double> v1, std::vector<double> v2);
double squaredEuclideanNorm(std::vector<double> v1);

void loadBinMatrix(Matrix &points, std::string filename);
void loadTxtMatrix(Matrix &points, std::string filename);
void loadTxtVector(std::vector<double> &vec, std::string filename);
void writeMatrix(const Matrix& data, std::string filename);
void writeVector(std::vector<double> data, std::string filename);

void print_vec(std::vector<double> vec);
void print_matrix(Matrix mat);

Matrix randSample(Matrix data, int num);
Matrix addZeroLoc(const Matrix& mat);
Matrix subMatrix(const Matrix& mat,int rowi,int coli,int rowf,int colf);
std::vector<double> subVector(std::vector<double> vec,int first, int end);
Matrix rescale(const Matrix& m, double min, double max);

#endif

