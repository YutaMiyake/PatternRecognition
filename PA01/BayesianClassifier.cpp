#include "BayesianClassifier.h"

BayesianClassifier::BayesianClassifier(){}
BayesianClassifier::~BayesianClassifier(){}

// Linear discriminant *****************************************************
LinearDiscriminant::LinearDiscriminant(std::vector<double> mean1, std::vector<double> mean2, double var, double prior1, double prior2){
  this->dim = mean1.size();
  this->mean1 = mean1;
  this->mean2 = mean2;
  this->var = var;
  this->prior1 = prior1;
  this->prior2 = prior2;
}

int LinearDiscriminant::predict(std::vector<double> x) const{
  double g1, g2;
  std::vector<double> diff1(this->dim), diff2(this->dim);
  double norm1, norm2;

  // x - mu
  for(int ith = 0; ith < this->dim; ith++){
    diff1[ith] = x[ith] - this->mean1[ith];
    diff2[ith] = x[ith] - this->mean2[ith];
  }

  // ||x - mu||^2
  norm1 = norm2 = 0.0;
  for(int ith = 0; ith < this->dim; ith++){
    norm1 += diff1[ith]*diff1[ith];
    norm2 += diff2[ith]*diff2[ith];
  }

  // p(w1) == p(w2)
  if(this->prior1 == this->prior2){
    g1 = -norm1;
    g2 = -norm2;
  }

  // p(w1) != p(w2)
  else{
    g1 = -norm1/(2*this->var*this->var) + log(this->prior1);
    g2 = -norm2/(2*this->var*this->var) + log(this->prior2);
  }

  if(g1 > g2){
    return 0;
  }
  return 1;
}

// Quadratic discriminant *****************************************************
QuadraticDiscriminant::QuadraticDiscriminant(std::vector<double> mean1, std::vector<double> mean2, Matrix cov1, Matrix cov2, double prior1, double prior2){
  this->dim = mean1.size();
  this->mean1 = mean1;
  this->mean2 = mean2;
  this->cov1 = cov1;
  this->cov2 = cov2;
  this->prior1 = prior1;
  this->prior2 = prior2;
  this->init();
}

void QuadraticDiscriminant::init(){
  // for now assume 2d distribution
  assert(this->mean1.size() == 2);
  assert(this->mean2.size() == 2);
  assert(this->cov1.size() == 2);
  assert(this->cov2.size() == 2);

  double cov1_inv[2][2], cov2_inv[2][2];
  double det1, det2;

  // determinant
  det1 = this->cov1[0][0]*this->cov1[1][1] - this->cov1[0][1]*this->cov1[1][0];
  det2 = this->cov2[0][0]*this->cov2[1][1] - this->cov2[0][1]*this->cov2[1][0];

  // calc inverse
  cov1_inv[0][0] = this->cov1[1][1]/det1;
  cov1_inv[0][1] = this->cov1[0][1]/det1;
  cov1_inv[1][0] = this->cov1[1][0]/det1;
  cov1_inv[1][1] = this->cov1[0][0]/det1;

  cov2_inv[0][0] = this->cov2[1][1]/det2;
  cov2_inv[0][1] = this->cov2[0][1]/det2;
  cov2_inv[1][0] = this->cov2[1][0]/det2;
  cov2_inv[1][1] = this->cov2[0][0]/det2;

  // W
  this->W1 = Matrix(2,std::vector<double>(2));
  this->W2 = Matrix(2,std::vector<double>(2));

  for(int row = 0; row < 2; row++){
    for(int col = 0; col < 2; col++){
      this->W1[row][col] = -0.5*cov1_inv[row][col];
      this->W2[row][col] = -0.5*cov2_inv[row][col];
    }
  }

  // w
  this->w1 = std::vector<double>(2,0.0);
  this->w2 = std::vector<double>(2,0.0);
  for(int row = 0; row < 2; row++){
    for(int col = 0; col < 2; col++){
      this->w1[row] += cov1_inv[row][col]*this->mean1[col];
      this->w2[row] += cov2_inv[row][col]*this->mean2[col];
    }
  }

  // w0
  double part11 = 0.0, part12 = 0.0;
  double part21 = 0.0, part22 = 0.0;

  std::vector<double> temp1(2,0.0), temp2(2,0.0);

  for(int col = 0; col < 2; col++){
    for(int row = 0; row < 2; row++){
      temp1[col] += cov1_inv[row][col]*this->mean1[row];
      temp2[col] += cov2_inv[row][col]*this->mean2[row];
    }
  }
  for(int row = 0; row < 2; row++){
      part11 += temp1[row]*this->mean1[row];
      part21 += temp2[row]*this->mean2[row];
  }
  part11 *= -0.5;
  part21 *= -0.5;

  det1 = cov1_inv[0][0]*cov1_inv[1][1] - cov1_inv[0][1]*cov1_inv[1][0];
  det2 = cov2_inv[0][0]*cov2_inv[1][1] - cov2_inv[0][1]*cov2_inv[1][0];

  part12 = -0.5*log(det1)+log(this->prior1);
  part22 = -0.5*log(det2)+log(this->prior2);
  this->w10 = part11 + part12;
  this->w20 = part21 + part22;

  assert(this->w1.size() == 2);
  assert(this->w2.size() == 2);
  assert(this->W1.size() == 2);
  assert(this->W2.size() == 2);
}

int QuadraticDiscriminant::predict(std::vector<double> x) const{
  std::vector<double> temp1(2,0.0), temp2(2,0.0);
  for(int col = 0; col < 2; col++){
    for(int row = 0; row < 2; row++){
      temp1[col] += this->W1[row][col]*x[row];
      temp2[col] += this->W2[row][col]*x[row];
    }
  }
  double dot11 = 0.0, dot21 = 0.0;
  for(int row = 0; row < 2; row++){
    dot11 += temp1[row]*x[row];
    dot21 += temp2[row]*x[row];
  }

  double dot12 = 0.0, dot22 = 0.0;
  for(int row = 0; row < 2; row++){
    dot12 += this->w1[row]*x[row];
    dot22 += this->w2[row]*x[row];
  }

  double g1, g2;
  g1 = dot11 + dot12 + this->w10;
  g2 = dot21 + dot22 + this->w20;

  if(g1 > g2){
    return 0;
  }
  return 1;
}