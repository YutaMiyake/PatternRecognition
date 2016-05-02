#include <iostream>
#include <fstream>
#include <sstream>
#include "libsvm-3.21/svm.h"
#include "Debugger.h"
#include "Matrix.h"
#include "ML.h"
#include "BayesianClassifier.h"

void ex1(int kernel_type, int degree, double cval, double gamma);
void ex2();

int main(){

  /*// polynomial svm
  int kernel_type = 1;
  for(int degree = 1; degree <=3; degree++){
    for(double cval = 1; cval <= 1000; cval*=10){
      ex1(kernel_type,degree,cval,1);
    }
  }

  // RBF svm
  kernel_type = 2;
  for(double gamma = 0.1; gamma >= 0.001; gamma/=10){
    for(double cval = 1; cval <= 1000; cval*=10){
      ex1(kernel_type,0,cval,gamma);
    }
  }
  */

  ex2();
  return 0;
}
void ex1(int kernel_type, int degree, double cval, double gamma){
  int train_size = 134;
  int valid_size = 133;
  int test_size = 133;
  int eigen_size = 30;

  // parameter initialization
  svm_parameter param;
  param.kernel_type = kernel_type;
  param.degree = degree;
  param.C = cval;
  param.svm_type = 0;
  param.gamma = gamma;
  param.coef0 = 0.0;
  param.nu = 0.5;
  param.cache_size = 100;
  param.eps = 1e-3;
  param.p = 0.1;
  param.shrinking = 1;
  param.probability = 0;
  param.nr_weight = 0;
  param.weight_label = NULL;
  param.weight = NULL;

  // read training data
  svm_node **node = new svm_node*[train_size];
  for(int data_idx = 0; data_idx < train_size; data_idx++){
    node[data_idx] = new svm_node[eigen_size + 1];
  }

  std::ifstream fin;
  std::string dummy;
  fin.clear();
  fin.open("genderdata/48_60/trPCA_03.txt");

  for(int eigen = 0; eigen < eigen_size; eigen++){
    for(int data_idx = 0; data_idx < train_size; data_idx++){
      node[data_idx][eigen].index = eigen + 1;
      fin >> node[data_idx][eigen].value;
    }
  }
  for(int data_idx = 0; data_idx < train_size; data_idx++){
    // Make delimiter node
    node[data_idx][eigen_size].index = -1;
    node[data_idx][eigen_size].value = 0.0;
  }
  fin.close();

  // read labels
  double* labels = new double[train_size];
  fin.clear();
  fin.open("genderdata/48_60/TtrPCA_03.txt");

  for(int label_idx = 0; label_idx < train_size; label_idx++){
    fin >> labels[label_idx];
  }
  fin.close();

  // set up svm problem
  svm_problem prob;
  prob.x = node;
  prob.y = labels;
  prob.l = train_size;

  // create model
  svm_model *model = svm_train(&prob, &param);

  // read test data (test section)
  svm_node **test_node = new svm_node*[test_size];
  for(int data_idx = 0; data_idx < train_size; data_idx++){
    test_node[data_idx] = new svm_node[eigen_size + 1];
  }

  fin.clear();
  fin.open("genderdata/48_60/tsPCA_03.txt");

  for(int eigen = 0; eigen < eigen_size; eigen++){
    for(int data_idx = 0; data_idx < test_size; data_idx++){
      test_node[data_idx][eigen].index = eigen + 1; // libsvm has no 0th data
      fin >> test_node[data_idx][eigen].value;
    }
  }
  for(int data_idx = 0; data_idx < train_size; data_idx++){
    // Make delimiter node
    test_node[data_idx][eigen_size].index = -1;
    test_node[data_idx][eigen_size].value = 0.0;
  }
  fin.close();

  // read test labels
  double* test_labels = new double[train_size];
  fin.clear();
  fin.open("genderdata/48_60/TtsPCA_03.txt");

  for(int label_idx = 0; label_idx < train_size; label_idx++){
    fin >> test_labels[label_idx];
  }
  fin.close();

  // test test data
  double predicted_label;
  double err_counter = 0;
  for(int data_idx = 0; data_idx < test_size; data_idx++){
    predicted_label = svm_predict(model, test_node[data_idx]);
    // std::cout << predicted_label << std::endl;
    if(predicted_label != test_labels[data_idx]){
      err_counter++;
    }
  }

  // read test data (validation section)
  svm_node **valid_node = new svm_node*[valid_size];
  for(int data_idx = 0; data_idx < train_size; data_idx++){
    valid_node[data_idx] = new svm_node[eigen_size + 1];
  }

  fin.clear();
  fin.open("genderdata/48_60/valPCA_03.txt");

  for(int eigen = 0; eigen < eigen_size; eigen++){
    int data_idx = 0;
    for(data_idx = 0; data_idx < valid_size; data_idx++)
    {
      valid_node[data_idx][eigen].index = eigen + 1;
      fin >> valid_node[data_idx][eigen].value;
    }
  }
  for(int data_idx = 0; data_idx < train_size; data_idx++){
    // Make delimiter node
    valid_node[data_idx][eigen_size].index = -1;
    valid_node[data_idx][eigen_size].value = 0.0;
  }
  fin.close();

  // read valid labels
  double* valid_labels = new double[valid_size];
  fin.clear();
  fin.open("genderdata/48_60/TvalPCA_03.txt");

  for(int label_idx = 0; label_idx < valid_size; label_idx++){
    fin >> valid_labels[label_idx];
  }
  fin.close();

  // test valid data
  for(int data_idx = 0; data_idx < valid_size; data_idx++){
    predicted_label = svm_predict(model, valid_node[data_idx]);
    if(predicted_label != valid_labels[data_idx]){
      err_counter++;
    }
  }

  Debugger logger("log_48_60",true);
  std::stringstream ss;
  if(kernel_type == POLY){
    ss << "Poly: degree = " << degree << " C = " << cval;
  }
  else{
    ss << "RBF: gamma = " << gamma << " C = " << cval;
  }
  ss << " err count : " << err_counter;
  ss << " Accuracy : " << 1 - err_counter/(test_size+valid_size);
  logger.debug(ss.str());

}
void ex2(){
  int male_size = 69;
  int female_size = 65;
  int train_size = male_size + female_size;
  int eigen_size = 30;

  // read labels
  double* labels = new double[train_size];
  std::ifstream fin;
  std::string dummy;
  fin.clear();
  fin.open("genderdata/48_60/TtrPCA_03.txt");
  for(int label_idx = 0; label_idx < train_size; label_idx++){
    fin >> labels[label_idx];
  }
  fin.close();

  // read train data
  Matrix train_data(train_size,std::vector<double>(eigen_size));
  Matrix male(male_size,std::vector<double>(eigen_size));
  Matrix female(female_size,std::vector<double>(eigen_size));

  fin.clear();
  fin.open("genderdata/48_60/trPCA_03.txt");

  for(int eigen = 0; eigen < eigen_size; eigen++){
    for(int data_idx = 0; data_idx < train_size; data_idx++){
      fin >> train_data[data_idx][eigen];
    }
  }
  fin.close();

  std::vector<double> temp(30);
  int midx = 0;
  int fidx = 0;
  for(int data_idx = 0; data_idx < train_size; data_idx++){
    if(labels[data_idx] == 1){
      male[midx++] = train_data[data_idx];
    }else if(labels[data_idx] == 2){
      female[fidx++] = train_data[data_idx];
    }
  }
  fin.close();

  // calc mean and covariance
  std::cout << "Calculating mean ..." << std::endl;
  std::vector<double> mean_male = getSampleMean(male);
  std::vector<double> mean_female = getSampleMean(female);

  std::cout << "Calculating covariance ..." << std::endl;
  Matrix cov_male = getSampleVar(male,mean_male);
  Matrix cov_female = getSampleVar(female,mean_female);

  // model discriminant
  std::cout << "Setting up model ..." << std::endl;
  QuadraticDiscriminant classifier(mean_male, mean_female, cov_male, cov_female);

  // read test data
  int test_size = 133;
  Matrix test_data(test_size,std::vector<double>(eigen_size));

  fin.clear();
  fin.open("genderdata/48_60/tsPCA_03.txt");

  for(int eigen = 0; eigen < eigen_size; eigen++){
    for(int data_idx = 0; data_idx < train_size; data_idx++){
      fin >> test_data[data_idx][eigen];
    }
  }

  fin.close();

  // read test labels
  double* test_labels = new double[test_size];
  fin.clear();
  fin.open("genderdata/48_60/TtsPCA_03.txt");
  for(int label_idx = 0; label_idx < train_size; label_idx++){
    fin >> labels[label_idx];
  }
  fin.close();

  // testing
  std::cout << "Testing model ..." << std::endl;
  int predicted_label;
  double err_ctr = 0;
  for(int data_idx = 0; data_idx < test_size; data_idx++){
    predicted_label = classifier.predict(test_data[data_idx]) + 1;
    if(predicted_label != test_labels[data_idx]){
      if(test_labels[data_idx] == 1){
        err_ctr++;
      }
      else{
        err_ctr++;
      }
    }
  }

  std::cout << "Err rate" << err_ctr/test_size << std::endl;

}
