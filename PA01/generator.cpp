#include "Generator.h"
Generator::Generator(){}
void Generator::generateSample(double mean, double std, int sampleNum, std::string filename){
  double xPos, yPos;
  for(int sample = 0; sample < sampleNum; sample++){
    xPos = generateGaussianNoise(mean,std);
    yPos = generateGaussianNoise(mean,std);
    writePointToFile(filename, xPos, yPos);
  }
}
void Generator::generateSample(double mean1, double mean2, double std1, double std2, int sampleNum, std::string filename){
  double xPos, yPos;
  for(int sample = 0; sample < sampleNum; sample++){
    xPos = generateGaussianNoise(mean1,std1);
    yPos = generateGaussianNoise(mean2,std2);
    writePointToFile(filename, xPos, yPos);
  }
}

void Generator::writePointToFile(std::string filename, double xPos, double yPos){
  std::ofstream fout;
  fout.open(filename.c_str(), std::ofstream::app);
  fout << xPos << " " << yPos << std::endl;
  fout.close();
}

double Generator::generateGaussianNoise(double mean, double std)
 /* Gaussian Random Number Generator. N(mean, std^2).
 It uses box-muller transformation to simulate the Gaussian distribution
 from uniformly distributed numbers.
 */
{
  double x1, x2, w, y1;
  static double y2;
  static bool use_last = false;

  // setup time-based seed
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> uni_real_dist(0.0,1.0);

  if (use_last){ /* use value from previous call */
    y1 = y2;
    use_last = false;
  }
  else{
    do {
      x1 = 2.0 * uni_real_dist(generator) - 1.0;
      x2 = 2.0 * uni_real_dist(generator) - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = true;
  }

  return( mean + y1 * std );
}