#ifndef __BAYESIAN_CLASSIFIER__H__
#define __BAYESIAN_CLASSIFIER__H__

#include <vector>
#include <cmath>
#include <cassert>
#include <iostream>
#include "Matrix.h"

class BayesianClassifier{
public:
  BayesianClassifier();
  virtual ~BayesianClassifier();
  virtual int predict(std::vector<double> x) const =0;

  double kb(double beta, std::vector<double> mean1, std::vector<double> mean2, Matrix cov1, Matrix cov2) const;
  double ekb(double beta) const;
  double chernoff(double beta) const;
  double bhattacharyya() const;

protected:
  int dim;
  std::vector<double> mean1, mean2;
  Matrix cov1, cov2;
  double prior1, prior2;
};

class LinearDiscriminant: public BayesianClassifier{
public:
  LinearDiscriminant(std::vector<double> mean1, std::vector<double> mean2, double var = 1.0, double prior1 = 0.5, double prior2 = 0.5);
  void init();
  int predict(std::vector<double> x) const;

private:
  double var;
};

class QuadraticDiscriminant: public BayesianClassifier{
public:
  QuadraticDiscriminant(std::vector<double> mean1, std::vector<double> mean2, Matrix cov1, Matrix cov2, double prior1 = 0.5, double prior2 = 0.5);
  void init();
  void setPriors(double prior1, double prior2);
  int predict(std::vector<double> x) const;

private:
  Matrix W1, W2;
  double w10, w20;
  std::vector<double> w1, w2;
};

#endif