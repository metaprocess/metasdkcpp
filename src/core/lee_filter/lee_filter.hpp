#include <Eigen/Dense>
#include <complex>
#include <cmath>
#include <algorithm>
#include "MatrixUtils.h"

Eigen::MatrixXcd LeeFilter(const Eigen::MatrixXcd& input_image, int window_size, int ENL);

// MatrixCompareResult test_lee_filter();