#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <string>

double box_muller(double m, double s);
void writePointToFile(std::string filename, double x, double y);
void generateSample(double mean, double std, int sampleNum, std::string filename);
void generateSample(double mean1, double mean2, double std1, double std2, int sampleNum, std::string filename);

int main(int argc, char** argv){
  // default setting
  std::string filename1 = "sample_data/output1";
  std::string filename2 = "sample_data/output2";
  std::string filename3 = "sample_data/output3";

  // remove file if exists
  std::remove(filename1.c_str());
  std::remove(filename2.c_str());
  std::remove(filename3.c_str());

  // generate part1 class1 and class 2
  generateSample(1.0, sqrt(1.0), 10000, filename1);
  generateSample(4.0, sqrt(1.0), 10000, filename2);

  // generate part2 class 2
  generateSample(4.0, 4.0, sqrt(4.0), sqrt(16.0), 10000, filename3);

   // Example
   std::remove("Example1");
   std::remove("Example2");
   generateSample(3.0, 6.0, sqrt(0.5), sqrt(2.0), 10000, "sample_data/Example1");
   generateSample(3.0, -2.0, sqrt(2.0), sqrt(2.0), 10000, "sample_data/Example2");
}

void generateSample(double mean, double std, int sampleNum, std::string filename){
  double xPos, yPos;
  for(int sample = 0; sample < sampleNum; sample++){
    xPos = box_muller(mean,std);
    yPos = box_muller(mean,std);
    writePointToFile(filename, xPos, yPos);
  }
}
void generateSample(double mean1, double mean2, double std1, double std2, int sampleNum, std::string filename){
  double xPos, yPos;
  for(int sample = 0; sample < sampleNum; sample++){
    xPos = box_muller(mean1,std1);
    yPos = box_muller(mean2,std2);
    writePointToFile(filename, xPos, yPos);
  }
}

void writePointToFile(std::string filename, double xPos, double yPos){
  std::ofstream fout;
  fout.open(filename.c_str(), std::ofstream::app);
  fout << xPos << " " << yPos << std::endl;
  fout.close();
}

double box_muller(double m, double s)
 /* normal random variate generator */
{
  /* mean m, standard deviation s */
  double x1, x2, w, y1;
  static double y2;
  static bool use_last = false;

  // setup time-based seed
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  if (use_last){ /* use value from previous call */
    y1 = y2;
    use_last = false;
  }
  else{
    do {
      x1 = 2.0 * distribution(generator) - 1.0;
      x2 = 2.0 * distribution(generator) - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = true;
  }

  return( m + y1 * s );
}