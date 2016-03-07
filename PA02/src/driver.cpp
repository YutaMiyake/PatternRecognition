#include "Generator.h"
#include "BayesianClassifier.h"
#include "ML.h"
#include <iostream>
#include <fstream>
#include "util.h"
#include "Debugger.h"
#include "Matrix.h"
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cassert>
#include "tools/ReadImage.cpp"
#include "tools/ReadImageHeader.cpp"
#include "tools/WriteImage.cpp"
#include "tools/rgb.h"

// function prototypes *****************************************************
void calcError(const BayesianClassifier &classifier,
  const Matrix &points1, const Matrix &points2, std::string filename);
void makeColorMatrices(ImageType& img, ImageType& ref, Matrix &sk_cols,
  Matrix &nsk_cols, bool YCbCr);
void testSkinRecognition(const BayesianClassifier &classifier, ImageType &image, ImageType &ref, std::string out, bool YCbCr);
void testSkinRecogWithThreshold(const std::vector<double> &mean, const Matrix &cov, ImageType &image, std::string out);

void genFiles();
void prob1a();
void prob1b();
void prob2a();
void prob2b();
void prob3a();
void prob3b();

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

  // problems
  // prob1a();
  // prob1b();
  // prob2a();
  // prob2b();
     prob3a();
  // prob3b();
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

  // Prob1: generate samples for class1 and class 2
  generator.generateSample(1.0, sqrt(2.0), 10000, filename1);
  generator.generateSample(6.0, sqrt(2.0), 10000, filename2);

  // Prob2: generate samples for class 2
  generator.generateSample(6.0, 6.0, sqrt(4.0), sqrt(8.0), 10000, filename3);
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

void makeColorMatrices(ImageType& img, ImageType& ref, Matrix &sk_cols,
  Matrix &nsk_cols, bool YCbCr)
{
  int height1, width1, levels1;
  int height2, width2, levels2;
  img.getImageInfo(height1, width1, levels1);
  ref.getImageInfo(height2, width2, levels2);

  assert(height1 == height2);
  assert(width1 == width2);
  assert(levels1 == levels2);

  RGB val1, val2;
  std::vector<double> color(2);
  RGB black(0,0,0);

  for(int row = 0; row < height1; row++){
    for(int col = 0; col < width1; col++){
      img.getPixelVal(row, col, val1);
      ref.getPixelVal(row, col, val2);

      if(YCbCr == true){
        color[0] = -0.169*val1.r - 0.332*val1.g+ 0.500*val1.b;
        color[1] = 0.500*val1.r - 0.419*val1.g - 0.081*val1.b;
      }
      else{
        color[0] = val1.r/float(val1.r+val1.g+val1.b);
        color[1] = val1.g/float(val1.r+val1.g+val1.b);
      }

      if(val2 != black){
        sk_cols.push_back(color);
      }
      else{
        nsk_cols.push_back(color);
      }
    }
  }
}

void testSkinRecognition(const BayesianClassifier &classifier, ImageType &image, ImageType &ref, std::string out, bool YCbCr){

    int height, width, levels;
    image.getImageInfo(height,width,levels);
    ImageType outImg(height,width,levels);

    RGB val1, val2;
    int label;
    std::vector<double> color(2);
    int TP = 0, TN = 0, FN = 0, FP = 0;
    RGB white(255,255,255);
    RGB black(0,0,0);

    for(int row = 0; row < height; row++){
      for(int col = 0; col < width; col++){
        image.getPixelVal(row, col, val1);
        ref.getPixelVal(row, col, val2);

        if(YCbCr == true){
          color[0] = -0.169*val1.r - 0.332*val1.g+ 0.500*val1.b;
          color[1] = 0.500*val1.r - 0.419*val1.g - 0.081*val1.b;
        }
        else{
          color[0] = val1.r/float(val1.r+val1.g+val1.b);
          color[1] = val1.g/float(val1.r+val1.g+val1.b);
        }

        label = classifier.predict(color);

        if(label == 0){
          outImg.setPixelVal(row, col, white);
          if(val2 != black){ TP++; } else{ FP++; }
        }
        else{
          outImg.setPixelVal(row, col, black);
          if(val2 == black){ TN++; } else{ FN++; }
        }
      }
    }  // end outer for loop

    std::cout << std::endl
              << "TP: " << TP << std::endl
              << "TN: " << TN << std::endl
              << "FP: " << FP << std::endl
              << "FN: " << FN << std::endl;

    /*std::stringstream ss;
    ss << FP << " " << FN;
    Debugger debugger("Data_Prog2/errors3a.txt",true);
    debugger.debug(ss.str());
    */

    writeImage(out.c_str(), outImg);
}
void testSkinRecogWithThreshold(const std::vector<double> &mean, const Matrix &cov, ImageType &image, std::string out){

    RGB white(255,255,255);
    RGB black(0,0,0);

    int height, width, levels;
    image.getImageInfo(height,width,levels);

    RGB val;
    std::vector<double> pc(2); // pure color

    double thR, thG;
    for(int row = 0; row < height; row++){
     for(int col = 0; col < width; col++){
       image.getPixelVal(row, col, val);

       pc[0] = val.r/float(val.r+val.g+val.b);
       pc[1] = val.g/float(val.r+val.g+val.b);

       thR = exp(-(cov[0][0] * pow((pc[0] - mean[0]),2) +  cov[0][1] * (pc[0]- mean[0])));
       thG = exp(-(cov[1][0] * (pc[1] - mean[1]) + cov[1][1] * pow((pc[1] - mean[1]),2)));

       if((thR >= .9 && thG >= 1.0 && thG < 1.2)
         || (thR <= .8 && thR >= .7 && thG > 1.1)){
         image.setPixelVal(row, col, white);
       }
       else{
         image.setPixelVal(row, col, black);
       }
      }
     } // end outer for loop

    writeImage(out.c_str(), image);
}

// problems ******************************************************************
void prob1a(){
  std::cout << "Problem 1 Part A ======== " << std::endl;

  Matrix points1(10000,std::vector<double>(2));
  Matrix points2(10000,std::vector<double>(2));
  loadFile(points1, "sample_data/output1");
  loadFile(points2, "sample_data/output2");

  std::vector<double> mean1 = getSampleMean(points1);
  std::vector<double> mean2 = getSampleMean(points2);

  std::cout << "sample mean1 = ";
  print_vec(mean1);
  std::cout << "sample mean2 = ";
  print_vec(mean2);

  Matrix cov1 = getSampleVar(points1, mean1);
  Matrix cov2 = getSampleVar(points2, mean2);

  std::cout << "sample cov1 = ";
  print_matrix(cov1);
  std::cout << "sample cov2 = ";
  print_matrix(cov2);

  std::cout << "With estimated params: " << std::endl;
  QuadraticDiscriminant classifier1(mean1, mean2, cov1, cov2);
  calcError(classifier1,points1,points2, "sample_data/labels1");

  std::cout << "Given params: " << std::endl;
  mean1 = {1.0,1.0};
  mean2 = {6.0,6.0};
  double var = 2.0;
  LinearDiscriminant classifier2(mean1, mean2, var);
  calcError(classifier2,points1,points2, "sample_data/labels1");
  std::cout << std::endl;
}
void prob1b(){
  std::cout << "Problem 1 Part B ======== " << std::endl;

  Matrix data1(10000,std::vector<double>(2));
  Matrix data2(10000,std::vector<double>(2));
  loadFile(data1, "sample_data/output1");
  loadFile(data2, "sample_data/output2");

  Matrix points1 = randSample(data1, 1000);
  Matrix points2 = randSample(data2, 1000);

  std::vector<double> mean1 = getSampleMean(points1);
  std::vector<double> mean2 = getSampleMean(points2);

  std::cout << "sample mean1 = ";
  print_vec(mean1);
  std::cout << "sample mean2 = ";
  print_vec(mean2);

  Matrix cov1 = getSampleVar(points1, mean1);
  Matrix cov2 = getSampleVar(points2, mean2);

  std::cout << "sample cov1 = ";
  print_matrix(cov1);
  std::cout << "sample cov2 = ";
  print_matrix(cov2);

  std::cout << "With estimated params: " << std::endl;
  QuadraticDiscriminant classifier1(mean1, mean2, cov1, cov2);
  calcError(classifier1,data1,data2, "sample_data/labels1");
  std::cout << std::endl;
}
void prob2a(){
  std::cout << "Problem 2 Part A ======== " << std::endl;

  Matrix points1(10000,std::vector<double>(2));
  Matrix points2(10000,std::vector<double>(2));
  loadFile(points1, "sample_data/output1");
  loadFile(points2, "sample_data/output3");

  std::vector<double> mean1 = getSampleMean(points1);
  std::vector<double> mean2 = getSampleMean(points2);

  std::cout << "sample mean1 = ";
  print_vec(mean1);
  std::cout << "sample mean2 = ";
  print_vec(mean2);

  Matrix cov1 = getSampleVar(points1, mean1);
  Matrix cov2 = getSampleVar(points2, mean2);

  std::cout << "sample cov1 = ";
  print_matrix(cov1);
  std::cout << "sample cov2 = ";
  print_matrix(cov2);

  std::cout << "With estimated params: " << std::endl;
  QuadraticDiscriminant classifier1(mean1, mean2, cov1, cov2);
  calcError(classifier1,points1,points2, "sample_data/labels1");

  std::cout << "Given params: " << std::endl;
  mean1 = {1.0,1.0};
  mean2 = {6.0,6.0};
  cov1 = {{2.0,0},{0,2.0}};
  cov2 = {{4.0,0},{0,8.0}};

  QuadraticDiscriminant classifier2(mean1, mean2, cov1, cov2);
  calcError(classifier2,points1,points2, "sample_data/labels1");
  std::cout << std::endl;
}
void prob2b(){
  std::cout << "Problem 2 Part B ======== " << std::endl;

  Matrix data1(10000,std::vector<double>(2));
  Matrix data2(10000,std::vector<double>(2));
  loadFile(data1, "sample_data/output1");
  loadFile(data2, "sample_data/output3");

  Matrix points1 = randSample(data1, 1000);
  Matrix points2 = randSample(data2, 1000);

  std::vector<double> mean1 = getSampleMean(points1);
  std::vector<double> mean2 = getSampleMean(points2);

  std::cout << "sample mean1 = ";
  print_vec(mean1);
  std::cout << "sample mean2 = ";
  print_vec(mean2);

  Matrix cov1 = getSampleVar(points1, mean1);
  Matrix cov2 = getSampleVar(points2, mean2);

  std::cout << "sample cov1 = ";
  print_matrix(cov1);
  std::cout << "sample cov2 = ";
  print_matrix(cov2);

  std::cout << "With estimated params: " << std::endl;
  QuadraticDiscriminant classifier1(mean1, mean2, cov1, cov2);
  calcError(classifier1,data1,data2, "sample_data/labels1");
  std::cout << std::endl;
}

void prob3a(){
  // filenames
  std::string train1 = "Data_Prog2/Training_1.ppm";
  std::string ref1 = "Data_Prog2/ref1.ppm";

  // variable declarations
  int M, N, Q;
  bool type;

  // make image objects
  readImageHeader(train1.c_str(), N, M, Q, type);
  ImageType image1(N, M, Q);
  ImageType refimage1(N, M, Q);
  readImage(train1.c_str(),image1);
  readImage(ref1.c_str(),refimage1);

  // make skin colors
  Matrix skin_colors, non_skin_colors;

  makeColorMatrices(image1, refimage1, skin_colors, non_skin_colors, false);

  // estimate parameters for skin-color class
  std::vector<double> mean1 = getSampleMean(skin_colors);
  std::cout << "sample mean1 = ";
  print_vec(mean1);
  Matrix cov1 = getSampleVar(skin_colors, mean1);
  std::cout << "sample cov1 = ";
  print_matrix(cov1);

  // set up parameters for non-skin-color class
  std::vector<double> mean2  = getSampleMean(non_skin_colors);
  Matrix cov2 =  getSampleVar(non_skin_colors, mean2);
  std::cout << "sample mean2 = ";
  print_vec(mean2);
  std::cout << "sample cov2 = ";
  print_matrix(cov2);

  // make classifier
  QuadraticDiscriminant classifier(mean1, mean2, cov1, cov2, 0.54, 0.46);

  // test skin classifications
  std::string out1 = "Data_Prog2/out1a.ppm";
  testSkinRecognition(classifier, image1, refimage1, out1, false);

  std::string train3 = "Data_Prog2/Training_3.ppm";
  std::string ref3 = "Data_Prog2/ref3.ppm";
  readImageHeader(train3.c_str(), N, M, Q, type);
  ImageType image3(N, M, Q);
  ImageType refimage3(N, M, Q);
  readImage(train3.c_str(),image3);
  readImage(ref3.c_str(),refimage3);
  std::string out3 = "Data_Prog2/out3a.ppm";
  testSkinRecognition(classifier, image3, refimage3, out3, false);

  std::string train6 = "Data_Prog2/Training_6.ppm";
  std::string ref6 = "Data_Prog2/ref6.ppm";
  readImageHeader(train6.c_str(), N, M, Q, type);
  ImageType image6(N, M, Q);
  ImageType refimage6(N, M, Q);
  readImage(train6.c_str(),image6);
  readImage(ref6.c_str(),refimage6);
  std::string out6 = "Data_Prog2/out6a.ppm";
  testSkinRecognition(classifier, image6, refimage6, out6, false);

}

void prob3b(){
  // filenames
  std::string train1 = "Data_Prog2/Training_1.ppm";
  std::string ref1 = "Data_Prog2/ref1.ppm";

  // variable declarations
  int M, N, Q;
  bool type;

  // make image objects
  readImageHeader(train1.c_str(), N, M, Q, type);
  ImageType image1(N, M, Q);
  ImageType refimage1(N, M, Q);
  readImage(train1.c_str(),image1);
  readImage(ref1.c_str(),refimage1);

  // make skin colors
  Matrix skin_colors, non_skin_colors;

  makeColorMatrices(image1, refimage1, skin_colors, non_skin_colors, true);

  // estimate parameters for skin-color class
  std::vector<double> mean1 = getSampleMean(skin_colors);
  std::cout << "sample mean1 = ";
  print_vec(mean1);
  Matrix cov1 = getSampleVar(skin_colors, mean1);
  std::cout << "sample cov1 = ";
  print_matrix(cov1);

  // set up parameters for non-skin-color class
  std::vector<double> mean2  = getSampleMean(non_skin_colors);
  Matrix cov2 =  getSampleVar(non_skin_colors, mean2);
  std::cout << "sample mean2 = ";
  print_vec(mean2);
  std::cout << "sample cov2 = ";
  print_matrix(cov2);

  // make classifier
  QuadraticDiscriminant classifier(mean1, mean2, cov1, cov2, 0.08, 0.92);

  std::string out1 = "Data_Prog2/out1b.ppm";
  testSkinRecognition(classifier, image1, refimage1, out1, true);

  std::string train3 = "Data_Prog2/Training_3.ppm";
  std::string ref3 = "Data_Prog2/ref3.ppm";
  readImageHeader(train3.c_str(), N, M, Q, type);
  ImageType image3(N, M, Q);
  ImageType refimage3(N, M, Q);
  readImage(train3.c_str(),image3);
  readImage(ref3.c_str(),refimage3);
  std::string out3 = "Data_Prog2/out3b.ppm";
  testSkinRecognition(classifier, image3, refimage3, out3, true);

  std::string train6 = "Data_Prog2/Training_6.ppm";
  std::string ref6 = "Data_Prog2/ref6.ppm";
  readImageHeader(train6.c_str(), N, M, Q, type);
  ImageType image6(N, M, Q);
  ImageType refimage6(N, M, Q);
  readImage(train6.c_str(),image6);
  readImage(ref6.c_str(),refimage6);
  std::string out6 = "Data_Prog2/out6b.ppm";
  testSkinRecognition(classifier, image6, refimage6, out6, true);
}