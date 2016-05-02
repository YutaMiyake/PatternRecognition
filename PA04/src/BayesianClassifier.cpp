#include "BayesianClassifier.h"

BayesianClassifier::BayesianClassifier(){}
BayesianClassifier::~BayesianClassifier(){}

// Error bounds **********************************************************
double BayesianClassifier::kb(double beta, std::vector<double> mean1, std::vector<double> mean2, Matrix cov1, Matrix cov2) const{

    std::vector<double> diff(this->dim);
    for(int ith = 0; ith < this->dim; ith++){
      diff[ith] = mean1[ith] - mean2[ith];
    }

    Matrix covline = cov1*(1-beta) + cov2*beta;

    double part1, part2;
    part1 =  beta*(1-beta)/2.0*dot_vecs(dot_mat_vec(inverse(covline),diff),diff);
    part2 = 0.5*(log(determinant(covline)) - log(pow(determinant(cov1),1-beta)*pow(determinant(cov2),beta)));

    double kb;
    kb = part1 + part2;
    return kb;
}
double BayesianClassifier::ekb(double beta) const{
  return exp(-kb(beta, this->mean1, this->mean2, this->cov1, this->cov2));
}

double BayesianClassifier::chernoff(double beta) const{
  return pow(this->prior1,beta)*pow(this->prior2,1-beta)*ekb(beta);
}

double BayesianClassifier::bhattacharyya() const{
  return sqrt(this->prior1*this->prior2)*ekb(0.5);
}

// Linear discriminant *****************************************************
LinearDiscriminant::LinearDiscriminant(std::vector<double> mean1, std::vector<double> mean2, double var, double prior1, double prior2){
  this->dim = mean1.size();
  this->mean1 = mean1;
  this->mean2 = mean2;
  this->var = var;
  this->prior1 = prior1;
  this->prior2 = prior2;
  this->init();
}

void LinearDiscriminant::init(){
  this->cov1 = Matrix(this->dim,std::vector<double>(this->dim));
  this->cov2 = Matrix(this->dim,std::vector<double>(this->dim));
  cov1[0] = {this->var,0};
  cov1[1] = {0,this->var};
  cov2[0] = {this->var,0};
  cov2[1] = {0,this->var};
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
  assert(this->mean1.size() == this->mean2.size());
  assert(this->cov1.size() == this->cov2.size());
  assert(this->cov1[0].size() == this->cov2[0].size());

  // calc inverse
  Matrix cov1_inv = inverse(this->cov1);
  Matrix cov2_inv = inverse(this->cov2);

  // W
  this->W1 = cov1_inv*(-0.5);
  this->W2 = cov2_inv*(-0.5);

  // w
  this->w1 = dot_mat_vec(cov1_inv,this->mean1);
  this->w2 = dot_mat_vec(cov2_inv,this->mean2);

  // w0
  double part11 = 0.0, part12 = 0.0;
  double part21 = 0.0, part22 = 0.0;

  part11 = -0.5*dot_vecs(this->w1,this->mean1);
  part21 = -0.5*dot_vecs(this->w2,this->mean2);

  part12 = -0.5*log(determinant(cov1_inv))+log(this->prior1);
  part22 = -0.5*log(determinant(cov2_inv))+log(this->prior2);
  this->w10 = part11 + part12;
  this->w20 = part21 + part22;
}

int QuadraticDiscriminant::predict(std::vector<double> x) const{
  double g1, g2;
  g1 = dot_vecs(dot_mat_vec(this->W1,x),x) + dot_vecs(this->w1,x) + this->w10;
  g2 = dot_vecs(dot_mat_vec(this->W2,x),x) + dot_vecs(this->w2,x) + this->w20;

  if(g1 > g2){
    return 0;
  }
  return 1;
}
void QuadraticDiscriminant::setPriors(double prior1, double prior2){
  this->prior1 = prior1;
  this->prior2 = prior2;
  init();
}