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
std::vector<double> vec_sum(std::vector<double> v1, std::vector<double> v2){
  assert(v1.size() == v2.size());
  int dim = v1.size();
  std::vector<double> temp(dim);
  for(int row = 0; row < dim; row++){
    temp[row] = v1[row] + v2[row];
  }
  return temp;
}
std::vector<double> vec_diff(std::vector<double> v1, std::vector<double> v2){
  assert(v1.size() == v2.size());
  int dim = v1.size();
  std::vector<double> temp(dim);
  for(int row = 0; row < dim; row++){
    temp[row] = v1[row] - v2[row];
  }
  return temp;
}

void loadFile(Matrix &points, std::string filename){
  int row = 0;
  int col = 0;
  double dummy;
  bool stop = false;
  std::ifstream fin;

  // input file
  fin.open(filename.c_str(), std::ifstream::binary);
  while(fin.good(), !stop){
    //std::cout << row << " " << col << std::endl;
    if(row < points.size()){
      if(col < points[row].size()){
        fin >> dummy;
        //std::cout << dummy << " ";
        points[row][col++] = dummy;
      }
      else{
        //std::cout << std::endl;
        col = 0;
        row++;
      }
    }
    else{
      stop = true;
    }
  }
  fin.close();
}
void print_vec(std::vector<double> vec){
  int size = vec.size();
  std::cout <<"[";
  for(int dim = 0; dim < size; dim++){
    std::cout << vec[dim];
    if(dim+1 < size){
      std::cout << ", ";
    }
  }
  std::cout << "]" << std::endl;
}
void print_matrix(Matrix mat){
  assert(!mat.empty());

  int row_size = mat.size();
  int col_size = mat[0].size();

  std::cout <<"[";
  for(int row = 0; row < row_size; row++){
    std::cout <<"[";
    for(int col = 0; col < col_size; col++){
      std::cout << mat[row][col];
      if(col+1 < col_size){
        std::cout << ", ";
      }
    }
    std::cout << "]";
    if(row+1 < row_size){
      std::cout << ", ";
    }
  }
  std::cout << "]" << std::endl;
}

Matrix randSample(Matrix data, int num){
  assert(data.size() >= num);

  int size = data.size();
  std::vector<double> indices(size);
  for(int idx = 0; idx < size; idx++){
    indices[idx] = idx;
  }
  std::random_shuffle ( indices.begin(), indices.end() );

  Matrix ret(num, std::vector<double>(data[0].size()));
  for(int point = 0; point < num; point++){
    ret[point] = data[indices[point]];
  }

  return ret;
}