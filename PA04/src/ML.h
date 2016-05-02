#ifndef __ML__H__
#define __ML__H__

#include <vector>
#include <cassert>
#include <iostream>
#include "Matrix.h"

std::vector<double> getSampleMean(const Matrix& train_data);
Matrix getSampleVar(const Matrix& train_data, const std::vector<double> &mean);

#endif