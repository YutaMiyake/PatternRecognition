#include "Generator.h"
#include "BayesianClassifier.h"
#include <iostream>
#include <fstream>

void errorTest(const BayesianClassifier &classifier,
  const Matrix &points1, const Matrix &points2);

void writeNumberToFile(std::string filename, double number);

void loadFile(Matrix &points, std::string filename);
void genFiles();
void prob1a();
void prob1b();
void prob2a();
void prob2b();
void prob3();

int main(int argvc, char** argv){
  // generate files
  genFiles();

  // classify
  prob1a();
  prob1b();
  prob2a();
  prob2b();
  prob3();
}

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

void errorTest(const BayesianClassifier &classifier,
  const Matrix &points1, const Matrix &points2){
  int label, error1 = 0, error2 = 0;

  for(int point = 0; point < points1.size(); point++){
   label = classifier.predict(points1[point]);
   if(label != 0){
    error1++;
   }
  }

  for(int point = 0; point < points2.size(); point++){
   label = classifier.predict(points2[point]);
   if(label != 1){
    error2++;
   }
  }
  std::cout << "Error classification for class1: " << error1 << std::endl;
  std::cout << "Error classification for class2: " << error2 << std::endl;
  std::cout << "Total error classification     : " << error1 + error2 << std::endl;
}

void writeNumberToFile(std::string filename, double number){
  std::ofstream fout;
  fout.open(filename.c_str(), std::ofstream::app);
  fout << number << std::endl;
  fout.close();
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
  LinearDiscriminant classifier1(mean1, mean2);

  errorTest(classifier1,points1,points2);
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
  errorTest(classifier2,points1,points2);
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
  errorTest(classifier3,points1,points2);
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
  errorTest(classifier4,points1,points2);
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
  errorTest(classifier5,points1,points2);
}