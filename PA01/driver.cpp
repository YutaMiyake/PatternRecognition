#include "Generator.h"
#include "BayesianClassifier.h"
#include <iostream>
#include <fstream>
#include "util.h"
#include "Debugger.h"
#include <sstream>
#include <iomanip>
#include <cstdlib>

// function prototypes *****************************************************
void calcError(const BayesianClassifier &classifier,
  const Matrix &points1, const Matrix &points2, std::string filename);

void loadFile(Matrix &points, std::string filename);
void genFiles();
void prob1a();
void prob1b();
void prob2a();
void prob2b();
void prob3();

// main function  *****************************************************
int main(int argvc, char** argv){

  // make directory
  const int dir_err = system("mkdir sample_data");
  if (dir_err == -1){
      printf("Error creating directory!n");
      exit(1);
  }

  // generate files
  genFiles();

  // classify
  prob1a();
  prob1b();
  prob2a();
  prob2b();
  prob3();
}

// supporting function implementation ***************************************
void genFiles(){
  Generator generator;

  // default filenames
  std::string filename1 = "sample_data/output1";
  std::string filename2 = "sample_data/output2";
  std::string filename3 = "sample_data/output3";

  // remove file if exists
  std::remove(filename1.c_str());
  std::remove(filename2.c_str());
  std::remove(filename3.c_str());

  // generate part1 class1 and class 2
  generator.generateSample(1.0, sqrt(2.0), 10000, filename1);
  generator.generateSample(6.0, sqrt(2.0), 10000, filename2);

  // generate part2 class 2
  generator.generateSample(6.0, 6.0, sqrt(4.0), sqrt(8.0), 10000, filename3);
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

void calcError(const BayesianClassifier &classifier,
  const Matrix &points1, const Matrix &points2, std::string filename){

  // for writing to files
  Debugger fileWriter(filename, true);
  std::remove(filename.c_str());
  std::stringstream ss;

  // calc errors
  int label, error1 = 0, error2 = 0;
  for(int point = 0; point < points1.size(); point++){
   label = classifier.predict(points1[point]);
   if(label != 0){
    error1++;
   }
   ss << label << "\n";
  }

  for(int point = 0; point < points2.size(); point++){
   label = classifier.predict(points2[point]);
   if(label != 1){
    error2++;
   }
   ss << label << "\n";
  }
  fileWriter.debug(ss.str());

  std::cout << "Error classification for class1: " << error1 << std::endl;
  std::cout << "Error classification for class2: " << error2 << std::endl;
  std::cout << "Total error classification     : " << error1 + error2
  << std::endl;
}

void calcErrorBounds(const BayesianClassifier &classifier, std::string filename){
  // for writing to files
  Debugger fileWriter(filename, true);
  std::remove(filename.c_str());
  std::stringstream ss;

  // minimize ekb
  double min_ekb, ekb, opt_beta, nth;
  std::vector<double> betas = util::linspace(0,1,1000);

  nth = 0;
  min_ekb = classifier.ekb(betas[nth]);

  for(nth = 1; nth < 1000; nth++){
    ekb = classifier.ekb(betas[nth]);
    if (min_ekb > ekb){
      min_ekb = ekb;
      opt_beta = betas[nth];
    }

    ss << betas[nth] << " " << ekb << "\n";
  }

  fileWriter.debug(ss.str());

  // outputs error bounds
  std::cout.precision(3);
  std::cout << "Opt-beta: " << opt_beta << std::endl
            << "Min-ekb : " << min_ekb << std::endl
            << "ekb(0.5) :" << classifier.ekb(0.5) << std::endl
            << "Chernoff bound: P(error) <= "
            << classifier.chernoff(opt_beta) << std::endl
            << "Bhattacharyya bound: P(error) <= "
            << classifier.bhattacharyya() << std::endl;
}
// problems ******************************************************************
void prob1a(){
  // minimum distance
  Matrix points1(10000,std::vector<double>(2));
  Matrix points2(10000,std::vector<double>(2));
  loadFile(points1, "sample_data/output1");
  loadFile(points2, "sample_data/output2");

  std::vector<double> mean1 = {1.0,1.0};
  std::vector<double> mean2 = {6.0,6.0};
  double var = 2.0;
  LinearDiscriminant classifier1(mean1, mean2, var);

  calcError(classifier1,points1,points2, "sample_data/labels1");
  calcErrorBounds(classifier1, "sample_data/error_bounds1");
  std::cout << std::endl;
}
void prob1b(){
  // linear discriminant
  Matrix points1(10000,std::vector<double>(2));
  Matrix points2(10000,std::vector<double>(2));
  loadFile(points1, "sample_data/output1");
  loadFile(points2, "sample_data/output2");

  std::vector<double> mean1 = {1.0,1.0};
  std::vector<double> mean2 = {6.0,6.0};
  double prior1 = 0.2;
  double prior2 = 0.8;
  double var = 2.0;

  LinearDiscriminant classifier2(mean1, mean2, var, prior1, prior2);
  calcError(classifier2,points1,points2, "sample_data/labels2");
  calcErrorBounds(classifier2, "sample_data/error_bounds2");
  std::cout << std::endl;
}
void prob2a(){
  // quadratic discriminant
  Matrix points1(10000,std::vector<double>(2));
  Matrix points2(10000,std::vector<double>(2));
  loadFile(points1, "sample_data/output1");
  loadFile(points2, "sample_data/output3");

  std::vector<double> mean1 = {1.0,1.0};
  std::vector<double> mean2 = {6.0,6.0};

  Matrix cov1(2,std::vector<double>(2));
  Matrix cov2(2,std::vector<double>(2));
  cov1[0] = {2.0,0};
  cov1[1] = {0,2.0};
  cov2[0] = {4.0,0};
  cov2[1] = {0,8.0};

  QuadraticDiscriminant classifier3(mean1, mean2, cov1, cov2);
  calcError(classifier3,points1,points2, "sample_data/labels3");
  calcErrorBounds(classifier3, "sample_data/error_bounds3");
  std::cout << std::endl;
}
void prob2b(){
  // quadratic discriminant
  Matrix points1(10000,std::vector<double>(2));
  Matrix points2(10000,std::vector<double>(2));
  loadFile(points1, "sample_data/output1");
  loadFile(points2, "sample_data/output3");

  std::vector<double> mean1 = {1.0,1.0};
  std::vector<double> mean2 = {6.0,6.0};

  Matrix cov1(2,std::vector<double>(2));
  Matrix cov2(2,std::vector<double>(2));
  cov1[0] = {2.0,0};
  cov1[1] = {0,2.0};
  cov2[0] = {4.0,0};
  cov2[1] = {0,8.0};

  double prior1 = 0.2;
  double prior2 = 0.8;

  QuadraticDiscriminant classifier4(mean1, mean2, cov1, cov2, prior1, prior2);
  calcError(classifier4,points1,points2, "sample_data/labels4");
  calcErrorBounds(classifier4, "sample_data/error_bounds4");
  std::cout << std::endl;
}
void prob3(){
  // quadratic discriminant
  Matrix points1(10000,std::vector<double>(2));
  Matrix points2(10000,std::vector<double>(2));
  loadFile(points1, "sample_data/output1");
  loadFile(points2, "sample_data/output3");

  std::vector<double> mean1 = {1.0,1.0};
  std::vector<double> mean2 = {6.0,6.0};

  LinearDiscriminant classifier5(mean1, mean2);
  calcError(classifier5,points1,points2, "sample_data/labels5");
  std::cout << std::endl;
}