#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <map>
#include "jacobi.h"
#include "ML.h"
#include "Matrix.h"

void expAI();
void expAII();
void expAV1();
void expAV2();
Matrix readData(bool train);
bool isMatched(std::vector<std::string> train_ids,
  std::vector<double> queryImage, std::string query_image_id, int numTopVec, int dim, std::string trainCoefFileName);
std::vector<std::string> readIds(std::string filename);

// constants
const int TRAIN_FILE_SIZE = 1204;
const int QUERY_FILE_SIZE = 1196;
const int WIDTH = 48;
const int HEIGHT = 60;
const int DIMENSION = WIDTH*HEIGHT;

int main(){
  // expAI();
  // expAII();
  expAV1();
  expAV2();
  return 0;
}

Matrix readData(bool train)
/*Read PGM images to matrix
  */
{
  std::stringstream ss;
  std::string dir = "Faces/fa_H/";
  std::string ext = ".pgm";
  std::ifstream fin;
  int file_size = TRAIN_FILE_SIZE;

  if(train == false){
    dir = "Faces/fb_H/";
    file_size = QUERY_FILE_SIZE;
  }

  Matrix ret(file_size, std::vector<double>(DIMENSION));

  // read image
  int val;
  char header[100];
  unsigned char * charImage;

  for(int fileNumber = 0; fileNumber < file_size; fileNumber++){
    ss << dir << fileNumber+1 << ext;

    fin.open(ss.str(), std::ios::in | std::ios::binary);
    if (!fin) {
      std::cout << "Can't read image: " << ss.str() << std::endl;
      exit(1);
    }

    fin.getline(header,100,'\n');
    if ( (header[0]!=80) ||    /* 'P' */
         (header[1]!=53) ) {   /* '5' */
         std::cout << "Image " << ss.str() << " is not PGM" << std::endl;
         exit(1);
    }

    charImage = (unsigned char *) new unsigned char [DIMENSION];
    fin.read( reinterpret_cast<char *>(charImage), (DIMENSION)*sizeof(unsigned char));

    for(int i=0; i<HEIGHT; i++){
       for(int j=0; j<WIDTH; j++){
         val = (int)charImage[i*WIDTH+j];
         ret[fileNumber][i*WIDTH+j] = val;
       }
    }

    ss.str("");
    fin.close();
  }
  return ret;
}

void writePGM(std::string filename, std::vector<double> img) {
    assert(img.size() == DIMENSION);

    std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary);

    fout << "P5" << WIDTH << " " << HEIGHT << " " << 255 << std::endl;

    unsigned char * charImage = (unsigned char *) new unsigned char [DIMENSION];
    int val;
    for(int i=0; i<HEIGHT; i++){
      for(int j=0; j<WIDTH; j++){
        val = img[i*WIDTH+j];
        charImage[i*WIDTH+j]=(unsigned char)val;
      }
    }

    fout.write( reinterpret_cast<char *>(charImage), (DIMENSION)*sizeof(unsigned char));
    fout.close();
}

void expAI(){
  /*
  Matrix dataMatrix = readData(true);

  // step 1 - calc mean
  std::vector<double> mean = getSampleMean(dataMatrix);
  writePGM("Faces/avgFace.pgm", mean);

  // step 2 - calc sample covariance
  // Matrix cov = getSampleVar(dataMatrix, mean);
  // writeMatrix(cov, "Faces/cov.txt");
  Matrix cov(DIMENSION,std::vector<double>(DIMENSION));
  loadTxtMatrix(cov, "Faces/cov_matrix.txt");

  // step 3 - calc eigenvectors and eigenvalues
  double** eigenvectors = NULL;
  double* eigenvalues = NULL;
  Matrix temp = addZeroLoc(cov);
  double** pcov = toPtr(temp);

  // allocate
  eigenvectors = new double*[DIMENSION+1];
  eigenvalues = new double[DIMENSION+1];
  for(int row = 0; row < DIMENSION+1; row++){
    eigenvectors[row] = new double[DIMENSION+1];
  }

  jacobi(pcov, DIMENSION, eigenvalues, eigenvectors);

  Matrix evecs = subMatrix(toMatrix(eigenvectors,DIMENSION+1,DIMENSION+1),1,1,DIMENSION+1,DIMENSION+1);
  std::vector<double> evals = subVector(toVector(eigenvalues,DIMENSION+1),1,DIMENSION+1);

  writeMatrix(evecs,"Faces/eigenvecs.txt");
  writeVector(evals,"Faces/eigenvals.txt");

  // deallocate
  for(int row = 0; row < DIMENSION+1; row++){
    delete[] eigenvectors[row];
  }
  delete[] eigenvalues;
  delete[] eigenvectors;

  // step 6 - normalize eigenvectors
  Matrix eigens(DIMENSION,std::vector<double>(DIMENSION));
  loadTxtMatrix(eigens, "Faces/eigenvecs.txt");
  eigens = transpose(eigens);

  for(int data = 0; data < TRAIN_FILE_SIZE + 1; data++){
    eigens[data] = normalize(eigens[data]);
  }
  writeMatrix(eigens,"Faces/normalized.txt");
  */

  Matrix normalized(DIMENSION,std::vector<double>(DIMENSION));
  loadTxtMatrix(normalized, "Faces/normalized.txt");

   // calc eigenfaces
   std::stringstream ss;
   std::string dir = "Faces/eigenfaces/";
   std::string ext = ".pgm";

   Matrix rescaled = rescale(normalized, 0, 255);

   // largest 10 eigenfaces
   for(int eigen = 0; eigen < 10; eigen++){
    ss << dir << eigen+1 << ext;
    writePGM(ss.str(),rescaled[eigen]);
    ss.str("");
   }

   // smallest 10 eigenfaces
   for(int eigen = DIMENSION-10; eigen < DIMENSION; eigen++){
    ss << dir << eigen+1 << ext;
    writePGM(ss.str(),rescaled[eigen]);
    ss.str("");
   }
}

void expAII(){
  Matrix trainDataMatrix = readData(true);

  std::vector<double> eigenvals(DIMENSION);
  loadTxtVector(eigenvals, "Faces/eigenvals.txt");

  // find k, preserving 80% of the info
  double info = 0.0;
  double sum1 = 0.0;
  double idx = 0;

  for(int row = 0; row < DIMENSION; row++){
    sum1 += eigenvals[row];
  }

  while(info/sum1 < 0.8){
    info += eigenvals[idx];
    idx++;
  }
  int dim = idx + 1;

  // calc eigen-coefficients
  /*
  Matrix eigens(DIMENSION,std::vector<double>(DIMENSION));
  loadTxtMatrix(eigens, "Faces/normalized.txt");

  std::vector<double> mean = getSampleMean(trainDataMatrix);

  Matrix trainCoeffs(TRAIN_FILE_SIZE,std::vector<double>(dim));
  for(int data = 0; data < TRAIN_FILE_SIZE; data++){
    trainDataMatrix[data] = trainDataMatrix[data] - mean; // subtract the avg face
    for(int nth = 0; nth < dim; nth++){
      trainCoeffs[data][nth] = dot_vecs(trainDataMatrix[data],eigens[nth]);
    }
  }
  writeMatrix(trainCoeffs,"Faces/trainCoeffs.txt");

  Matrix queryCoeffs(QUERY_FILE_SIZE,std::vector<double>(dim));
  Matrix queryDataMatrix = readData(false);
  for(int data = 0; data < QUERY_FILE_SIZE; data++){
    queryDataMatrix[data] = queryDataMatrix[data] - mean; // subtract the avg face
    for(int nth = 0; nth < dim; nth++){
      queryCoeffs[data][nth] = dot_vecs(queryDataMatrix[data],eigens[nth]);
    }
  }
  writeMatrix(queryCoeffs,"Faces/queryCoeffs.txt");
  */

  // read filenames
  std::vector<std::string> faH_ids(TRAIN_FILE_SIZE);
  faH_ids = readIds("Faces/fa_H_fileNames.txt");
  std::vector<std::string> fbH_ids(QUERY_FILE_SIZE);
  fbH_ids = readIds("Faces/fb_H_fileNames.txt");

  // testing
  Matrix queryCoeffs(QUERY_FILE_SIZE,std::vector<double>(dim));
  loadTxtMatrix(queryCoeffs, "Faces/queryCoeffs.txt");

  // Loop from N = 1 to 50
  for(int numTopVec = 1; numTopVec <= 50; numTopVec++){
    int count = 0;
    bool matched;
    for(int query = 0; query < QUERY_FILE_SIZE; query++){
      matched = isMatched(faH_ids, queryCoeffs[query], fbH_ids[query],numTopVec, dim, "Faces/trainCoeffs.txt");
      if(matched){
        count++;
        // std::cout << "Query " << query << " is correctly matched" << std::endl;
      }
    }
     std::cout << "N = " << numTopVec << " # of correctly matched query images = " << count
     << " performance = " << count/(double)QUERY_FILE_SIZE << std::endl;
  }
}

void expAV1(){
  Matrix trainDataMatrix = readData(true);

  std::vector<double> eigenvals(DIMENSION);
  loadTxtVector(eigenvals, "Faces/eigenvals.txt");

  // find k, preserving 90% of the info
  double info = 0.0;
  double sum1 = 0.0;
  double idx = 0;

  for(int row = 0; row < DIMENSION; row++){
    sum1 += eigenvals[row];
  }

  while(info/sum1 < 0.9){
    info += eigenvals[idx];
    idx++;
  }
  int dim = idx + 1;

  // calc eigen-coefficients
  /*
  Matrix eigens(DIMENSION,std::vector<double>(DIMENSION));
  loadTxtMatrix(eigens, "Faces/normalized.txt");

  std::vector<double> mean = getSampleMean(trainDataMatrix);

  Matrix trainCoeffs(TRAIN_FILE_SIZE,std::vector<double>(dim));
  for(int data = 0; data < TRAIN_FILE_SIZE; data++){
    trainDataMatrix[data] = trainDataMatrix[data] - mean; // subtract the avg face
    for(int nth = 0; nth < dim; nth++){
      trainCoeffs[data][nth] = dot_vecs(trainDataMatrix[data],eigens[nth]);
    }
  }
  writeMatrix(trainCoeffs,"Faces/trainCoeffs_90.txt");

  Matrix queryCoeffs(QUERY_FILE_SIZE,std::vector<double>(dim));
  Matrix queryDataMatrix = readData(false);
  for(int data = 0; data < QUERY_FILE_SIZE; data++){
    queryDataMatrix[data] = queryDataMatrix[data] - mean; // subtract the avg face
    for(int nth = 0; nth < dim; nth++){
      queryCoeffs[data][nth] = dot_vecs(queryDataMatrix[data],eigens[nth]);
    }
  }
  writeMatrix(queryCoeffs,"Faces/queryCoeffs_90.txt");
  */

  // read filenames
  std::vector<std::string> faH_ids(TRAIN_FILE_SIZE);
  faH_ids = readIds("Faces/fa_H_fileNames.txt");
  std::vector<std::string> fbH_ids(QUERY_FILE_SIZE);
  fbH_ids = readIds("Faces/fb_H_fileNames.txt");

  // testing
  Matrix queryCoeffs(QUERY_FILE_SIZE,std::vector<double>(dim));
  loadTxtMatrix(queryCoeffs, "Faces/queryCoeffs_90.txt");

  // Loop from N = 1 to 50
  for(int numTopVec = 1; numTopVec <= 50; numTopVec++){
    int count = 0;
    bool matched;
    for(int query = 0; query < QUERY_FILE_SIZE; query++){
      matched = isMatched(faH_ids, queryCoeffs[query], fbH_ids[query],numTopVec, dim, "Faces/trainCoeffs_90.txt");
      if(matched){
        count++;
        // std::cout << "Query " << query << " is correctly matched" << std::endl;
      }
    }
     std::cout << "N = " << numTopVec << " # of correctly matched query images = " << count
     << " performance = " << count/(double)QUERY_FILE_SIZE << std::endl;
  }

}

void expAV2(){
  Matrix trainDataMatrix = readData(true);

  std::vector<double> eigenvals(DIMENSION);
  loadTxtVector(eigenvals, "Faces/eigenvals.txt");

  // find k, preserving 95% of the info
  double info = 0.0;
  double sum1 = 0.0;
  double idx = 0;

  for(int row = 0; row < DIMENSION; row++){
    sum1 += eigenvals[row];
  }

  while(info/sum1 < 0.95){
    info += eigenvals[idx];
    idx++;
  }
  int dim = idx + 1;

  // calc eigen-coefficients
  /*
  Matrix eigens(DIMENSION,std::vector<double>(DIMENSION));
  loadTxtMatrix(eigens, "Faces/normalized.txt");

  std::vector<double> mean = getSampleMean(trainDataMatrix);

  Matrix trainCoeffs(TRAIN_FILE_SIZE,std::vector<double>(dim));
  for(int data = 0; data < TRAIN_FILE_SIZE; data++){
    trainDataMatrix[data] = trainDataMatrix[data] - mean; // subtract the avg face
    for(int nth = 0; nth < dim; nth++){
      trainCoeffs[data][nth] = dot_vecs(trainDataMatrix[data],eigens[nth]);
    }
  }
  writeMatrix(trainCoeffs,"Faces/trainCoeffs_95.txt");

  Matrix queryCoeffs(QUERY_FILE_SIZE,std::vector<double>(dim));
  Matrix queryDataMatrix = readData(false);
  for(int data = 0; data < QUERY_FILE_SIZE; data++){
    queryDataMatrix[data] = queryDataMatrix[data] - mean; // subtract the avg face
    for(int nth = 0; nth < dim; nth++){
      queryCoeffs[data][nth] = dot_vecs(queryDataMatrix[data],eigens[nth]);
    }
  }
  writeMatrix(queryCoeffs,"Faces/queryCoeffs_95.txt");
  */

  // read filenames
  std::vector<std::string> faH_ids(TRAIN_FILE_SIZE);
  faH_ids = readIds("Faces/fa_H_fileNames.txt");
  std::vector<std::string> fbH_ids(QUERY_FILE_SIZE);
  fbH_ids = readIds("Faces/fb_H_fileNames.txt");

  // testing
  Matrix queryCoeffs(QUERY_FILE_SIZE,std::vector<double>(dim));
  loadTxtMatrix(queryCoeffs, "Faces/queryCoeffs_95.txt");

  // Loop from N = 1 to 50
  for(int numTopVec = 1; numTopVec <= 50; numTopVec++){
    int count = 0;
    bool matched;
    for(int query = 0; query < QUERY_FILE_SIZE; query++){
      matched = isMatched(faH_ids, queryCoeffs[query], fbH_ids[query],numTopVec, dim, "Faces/trainCoeffs_95.txt");
      if(matched){
        count++;
        // std::cout << "Query " << query << " is correctly matched" << std::endl;
      }
    }
     std::cout << "N = " << numTopVec << " # of correctly matched query images = " << count
     << " performance = " << count/(double)QUERY_FILE_SIZE << std::endl;
  }

}

bool isMatched(std::vector<std::string> train_ids,
  std::vector<double> queryImage, std::string query_image_id, int numTopVec, int dim, std::string trainCoefFileName){
  // assume each distance is discrete
  Matrix trainCoeffs(TRAIN_FILE_SIZE,std::vector<double>(dim));
  loadTxtMatrix(trainCoeffs, trainCoefFileName);

  std::map<double,int> candidates;
  double dist;

  for(int train = 0; train < TRAIN_FILE_SIZE; train++){
    dist = squaredEuclideanNorm(trainCoeffs[train] - queryImage);
    candidates[dist] = train;
  }
  int ctr = 0;
  for (std::map<double,int>::iterator iter = candidates.begin();
    iter!=candidates.end() && ctr < numTopVec; ctr++, ++iter){
    if(train_ids[iter->second] == query_image_id){
      // std::cout << query_image_id << std::endl;
      return true;
    }
  }
  return false;
}

std::vector<std::string> readIds(std::string filename){
  std::fstream fin;
  fin.open(filename);
  std::vector<std::string> names;
  std::string name;
  while(fin.good()){
    fin >> name;
    names.push_back(name.substr(0,5));
  }
  fin.close();
  return names;
}

void expAV(){

}


