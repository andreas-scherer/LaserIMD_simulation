#pragma once
#include <Eigen/Eigen>
#include <complex>
using cfloat = std::complex<double>;
using Matrix4cd = Eigen::Matrix<cfloat, 4, 4>;
using Vector4cd = Eigen::Matrix<cfloat, 4, 1>;
using Matrix5cd = Eigen::Matrix<cfloat, 5, 5>;
using Matrix6cd = Eigen::Matrix<cfloat, 6, 6>;
using Vector6cd = Eigen::Matrix<cfloat, 6, 1>;
using Matrix8cd = Eigen::Matrix<cfloat, 8, 8>;
using Matrix9cd = Eigen::Matrix<cfloat, 9, 9>;
using Matrix12cd = Eigen::Matrix<cfloat, 12, 12>;
using Matrix18cd = Eigen::Matrix<cfloat, 18, 18>;
using Matrix24cd = Eigen::Matrix<cfloat, 24, 24>;
using Matrix36cd = Eigen::Matrix<cfloat, 36, 36>;
using Matrix64cd = Eigen::Matrix<cfloat, 64, 64>;
using Vector64cd = Eigen::Matrix<cfloat, 64, 1>;