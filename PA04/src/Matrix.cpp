#include "Matrix.h"

Matrix operator+(const Matrix& A, const Matrix& B)
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

Matrix operator*(const Matrix& A, double a)
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
Matrix operator/(const Matrix& A, double a){
  int nrows = A.size();
  int ncols = A[0].size();

  Matrix matrix(nrows, std::vector<double>(ncols));
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
    matrix[i][j] = A[i][j]/a;
    }
  }
  return matrix;
}
Matrix matrix_mul_matrix(const Matrix& A,const Matrix& B)
/* return AB */
{
  assert(!A.empty());
  assert(!B.empty());
  assert(A[0].size() == B.size());

  int cols = B[0].size();
  int rows = A.size();

  Matrix ret(rows, std::vector<double>(cols));
  Matrix BT = transpose(B);

  for(int row = 0; row < rows; row++){
    for(int col = 0; col < cols; col++){
      ret[row][col] = dot_vecs(A[row],BT[col]);
    }
  }
  return ret;
}
Matrix inverse(const Matrix& A){
  int nrows = A.size();
  int ncols = A[0].size();
  assert(nrows == ncols); // square matrix only

  std::cout << "Calculating inverse matrix (" << nrows << "x" << ncols << ") ..." <<std::endl;

  double det = 1.0/determinant(A);
  Matrix inv(nrows,std::vector<double>(ncols));
  Matrix cofactor(nrows-1,std::vector<double>(ncols-1));

  int cofactorCol,cofactorRow;

  for(int row = 0; row < ncols; row++){
    for(int col = 0; col < nrows; col++){
      cofactorRow = 0;
      for(int row2 = 0; row2 < nrows; row2++){
        if(row2 != row){
          cofactorCol = 0;
          for(int col2 = 0; col2 < ncols; col2++){
            if( col2 != col ){
              cofactor[cofactorCol][cofactorRow] = A[row2][col2];
              cofactorCol++;
            }
          }
          cofactorRow++;
        }
      }
      std::cout << '.';
      inv[col][row] = ((col+row)%2==1?-1.0:1.0)*det*determinant(cofactor);
    }
  }

  return inv;
}

Matrix matrix_vec_diff(const Matrix& A, std::vector<double> v){
  assert(A.size() >=1);
  assert(v.size() >=1);
  assert(A[0].size() == v.size());

  int rows = A.size();
  int cols = v.size();
  Matrix ret(rows,std::vector<double>(cols));
  for(int row = 0; row < rows; row++){
    for(int col = 0; col < cols; col++){
      ret[row][col] = A[row][col] - v[col];
    }
  }
  return ret;
}

Matrix transpose(const Matrix& A){
  assert(A.size() >= 1);
  assert(A[0].size() >= 1);

  int rows = A.size();
  int cols = A[0].size();
  Matrix transposed(cols, std::vector<double>(rows));
  for(int row = 0; row < rows; row++){
    for(int col = 0; col < cols; col++){
      transposed[col][row] = A[row][col];
    }
  }
  return transposed;
}

double determinant(const Matrix& A){
  int nrows = A.size();
  int ncols = A[0].size();
  assert(nrows == ncols); // only square matrix

  // base case
  if(nrows == 1 && ncols == 1){
    return A[0][0];
  }

  double det = 0;
  int cofactorCol, cofactorRow;
  int row, col;
  Matrix cofactor(nrows-1,std::vector<double>(ncols-1));

  for(col = 0; col < ncols; col++){
    // get next cofactor matrix
    cofactorRow = 0;
    for(row = 1; row < nrows; row++){
        cofactorCol = 0;
        for(int col2 = 0; col2 < ncols; col2++){
          if( col2 != col ){
            cofactor[cofactorCol][cofactorRow] = A[row][col2];
            cofactorCol++;
          }
        }
        cofactorRow++;
    }
    // recursive call
    det += (col%2==1?-1.0:1.0) * A[0][col] * determinant(cofactor);
  }

  return det;
}

Matrix toMatrix(double** A, int rows, int cols){
  assert(A != NULL);

  Matrix ret(rows,std::vector<double>(cols));

  for(int row = 0; row < rows; row++){
    for(int col = 0; col < cols; col++){
      ret[row][col] = A[row][col];
    }
  }
  return ret;
}

double** toPtr(const Matrix& A){
  assert(!A.empty() && !A[0].empty());

  int rows = A.size();
  int cols = A[0].size();

  double** ret = new double* [rows];

  for(int row = 0; row < rows; row++){
    ret[row] = new double [cols];
    for(int col = 0; col < cols; col++){
      ret[row][col] = A[row][col];
    }
  }
  return ret;
}

std::vector<double> dot_mat_vec(Matrix A, std::vector<double> v)
/* return Av */
{
  assert(A[0].size() == v.size());

  int rows = A.size();
  int cols = A[0].size();
  std::vector<double> dot(rows);

  for(int row = 0; row < rows; row++){
    for(int col = 0; col < cols; col++){
      dot[row] += A[row][col]*v[col];
    }
  }
  return dot;
}

double dot_vecs(std::vector<double> v,std::vector<double> u){
  assert(v.size() == u.size());
  double dot = 0.0;
  int dim = v.size();

  for(int row = 0; row < dim; row++){
    dot += v[row] * u[row];
  }
  return dot;
}

std::vector<double> operator+(std::vector<double> v1, std::vector<double> v2){
  assert(v1.size() == v2.size());
  int dim = v1.size();
  std::vector<double> temp(dim);
  for(int row = 0; row < dim; row++){
    temp[row] = v1[row] + v2[row];
  }
  return temp;
}

std::vector<double> operator-(std::vector<double> v1, std::vector<double> v2){
  assert(v1.size() == v2.size());
  int dim = v1.size();

  std::vector<double> temp(dim);
  for(int row = 0; row < dim; row++){
    temp[row] = v1[row] - v2[row];
  }
  return temp;
}
std::vector<double> toVector(double* A, int size){
  assert(A!=NULL);

  std::vector<double> ret(size);
  for(int row = 0; row < size; row++){
    ret[row] = A[row];
  }
  return ret;
}

std::vector<double> operator*(std::vector<double> v1, double a){
  int size = v1.size();
  std::vector<double> ret(size);

  for(int dim = 0; dim < size; dim++){
    ret[dim] = v1[dim] * a;
  }
  return ret;
}
double euclideanNorm(std::vector<double> v1){
  double sum = 0;
  int size = v1.size();
  for(int dim = 0; dim < size; dim++){
    sum += pow(v1[dim],2);
  }
  return sqrt(sum);
}
double euclideanDistance(std::vector<double> v1, std::vector<double> v2){
  assert(v1.size() == v2.size());
  std::vector<double> diff = v1-v2;
  double dist = euclideanNorm(diff);
  return dist;
}
double squaredEuclideanNorm(std::vector<double> v1){
  double sum = 0;
  int size = v1.size();
  for(int dim = 0; dim < size; dim++){
    sum += pow(v1[dim],2);
  }
  return sum;
}
std::vector<double> normalize(std::vector<double> v1){
  int size = v1.size();
  double length = euclideanNorm(v1);
  std::vector<double> ret(size);

  for(int dim = 0; dim < size; dim++){
   ret[dim] = v1[dim] / length;
  }
  return ret;
}

void loadBinMatrix(Matrix &points, std::string filename){
  int row = 0;
  int col = 0;
  double dummy;
  bool stop = false;
  std::ifstream fin;

  fin.open(filename.c_str(), std::ifstream::binary);

  if(!fin){
    std::cout << "Can't open " << filename << std::endl;
  }

  while(fin.good(), !stop){
    if(row < points.size()){
      if(col < points[row].size()){
        fin >> dummy;
        points[row][col++] = dummy;
      }
      else{
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
void loadTxtMatrix(Matrix &points, std::string filename){
  std::ifstream fin;
  fin.open(filename.c_str());
  if(!fin){
    std::cout << "Can't open " << filename << std::endl;
  }

  int rows, cols;
  double val;
  fin >> rows >> cols;
  for(int row = 0; row < rows; row++){
    for(int col = 0; col < cols; col++){
      fin >> val;
      points[row][col] = val;
    }
  }
  fin.close();
}
void loadTxtVector(std::vector<double> &vec, std::string filename){
  std::ifstream fin;
  fin.open(filename.c_str());
  if(!fin){
    std::cout << "Can't open " << filename << std::endl;
  }

  int rows;
  double val;
  fin >> rows;
  for(int row = 0; row < rows; row++){
    fin >> val;
    vec[row] = val;
  }

  fin.close();
}
void writeMatrix(const Matrix& data, std::string filename){
  std::ofstream fout;
  int rows = data.size();
  int cols = data[0].size();

  fout.open(filename);
  fout << rows << " " << cols << std::endl;
  for(int row = 0; row < rows; row++){
    for(int col = 0; col < cols; col++){
      fout << data[row][col];
      if(col+1 != cols){
        fout << " ";
      }
    }
    fout << std::endl;
  }
}
void writeVector(std::vector<double> data, std::string filename){
  std::ofstream fout;
  int rows = data.size();

  fout.open(filename);
  fout << rows << std::endl;
  for(int row = 0; row < rows; row++){
    fout << data[row];
    if(row+1 < rows){
      fout << std::endl;
    }
  }
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
Matrix addZeroLoc(const Matrix& mat) {
  int rows = mat.size()+1;
  int cols = mat[0].size()+1;

  Matrix ret(rows,std::vector<double>(cols));

  for(int row = 0; row < rows; row++){
    ret[row][0] = 0;
  }
  for(int col = 0; col < cols; col++){
    ret[0][col] = 0;
  }

  for(int row = 1; row < rows; row++){
    for(int col = 1; col < cols; col++){
      ret[row][col] = mat[row-1][col-1];
    }
  }
  return ret;
}
Matrix subMatrix(const Matrix& mat,int rowi,int coli,int rowf,int colf){
  assert(mat.size() >= 1);
  assert(mat[0].size() >= 1);
  assert(rowi >= 0 && coli >= 0);
  assert(rowf <= mat.size() && colf <= mat[0].size());

  Matrix ret(rowf-rowi, std::vector<double>(colf-coli));

  for(int row = rowi; row < rowf; row++){
    for(int col = coli; col < colf; col++){
      ret[row-rowi][col-coli] = mat[row][col];
    }
  }
  return ret;
}
std::vector<double> subVector(std::vector<double> vec, int first, int end){
  assert(vec.size() >= 1);
  assert(first >= 0 && end <= vec.size());
  assert(first <= end);
  std::vector<double> ret(end-first);

  for(int row = first; row < end; row++){
    ret[row-first] = vec[row];
  }
  return ret;
}
Matrix rescale(const Matrix& m, double min, double max)
/* rescale a value in the new range [max, min]*/
{
    int rows = m.size();
    int cols = m[0].size();
    double m_min = m[0][0];
    double m_max = m[0][0];
    for (int row = 0; row < rows; row++) {
        for (int col = 0; col < cols; col++) {
            if (m[row][col] < m_min) {
                m_min = m[row][col];
            }
            if (m[row][col] > m_max) {
                m_max = m[row][col];
            }
        }
    }
    double old_range = m_max - m_min;
    double new_range = max - min;

    // Create a new matrix with scaled elements.
    Matrix ret(rows, std::vector<double>(cols));
    for (int row = 0; row < rows; row++) {
      for (int col = 0; col < cols; col++) {
        ret[row][col] = ((m[row][col] - m_min)*new_range/old_range + min);
      }
    }
    return ret;
}