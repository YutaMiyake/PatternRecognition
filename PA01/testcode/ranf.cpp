#include <iostream>
#include <random>
int main(){
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  std::cout << distribution(generator) << std::endl;
}