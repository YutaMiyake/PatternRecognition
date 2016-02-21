#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <string>

class Generator{
public:
  Generator();
  double generateGaussianNoise(double m, double s);
  void writePointToFile(std::string filename, double x, double y);
  void generateSample(double mean, double std, int sampleNum, std::string filename);
  void generateSample(double mean1, double mean2, double std1, double std2, int sampleNum, std::string filename);
private:
};
