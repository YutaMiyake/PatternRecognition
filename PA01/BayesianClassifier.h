#include <vector>
#include <cmath>
#include <cassert>
#include <iostream>

typedef std::vector< std::vector<double> > Matrix;

class BayesianClassifier{
public:
  BayesianClassifier();
  virtual ~BayesianClassifier();
  virtual int predict(std::vector<double> x) const =0;
};

class LinearDiscriminant: public BayesianClassifier{
public:
  LinearDiscriminant(std::vector<double> mean1, std::vector<double> mean2, double var = 1.0, double prior1 = 0.5, double prior2 = 0.5);

  int predict(std::vector<double> x) const;

private:
  int dim;
  std::vector<double> mean1;
  std::vector<double> mean2;
  double var;
  double prior1;
  double prior2;
};

class QuadraticDiscriminant: public BayesianClassifier{
public:
  QuadraticDiscriminant(std::vector<double> mean1, std::vector<double> mean2, Matrix cov1, Matrix cov2, double prior1 = 0.5, double prior2 = 0.5);
  void init();
  int predict(std::vector<double> x) const;

private:
  int dim;
  std::vector<double> mean1, mean2;
  std::vector<double> w1, w2;
  Matrix cov1, cov2;
  double prior1, prior2;
  Matrix W1, W2;
  double w10, w20;
};