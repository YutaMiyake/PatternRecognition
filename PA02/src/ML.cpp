#include "ML.h"
#include <cmath>

std::vector<double> getSampleMean(const Matrix &train_data){
  assert(!train_data.empty());
  assert(!train_data[0].empty());

  int dim = train_data[0].size();
  int train_size = train_data.size();
  std::vector<double> mean(dim);

  for(int point = 0; point < train_size; point++){
    for(int nth = 0; nth < dim; nth++){
      if(!std::isnan(train_data[point][nth])){
        mean[nth] += train_data[point][nth] / train_size;
      }
    }
  }
  return mean;
}
Matrix getSampleVar(const Matrix &train_data, const std::vector<double>& mean){
  assert(!train_data.empty());
  assert(!train_data[0].empty());
  assert(!mean.empty());
  assert(train_data[0].size() == mean.size());

  int dim = train_data[0].size();
  int train_size = train_data.size();
  Matrix cov(dim, std::vector<double>(dim));

  std::vector<double> temp(dim);

  for(int point = 0; point < train_size; point++){
    temp = vec_diff(train_data[point], mean);
    for(int col = 0; col < dim; col++){
      for(int row = 0; row < dim; row++){
        if(!std::isnan(temp[row])){
          cov[row][col] +=  temp[row] * temp[col] / train_size;
        }
      }
    }
  }

  return cov;
}