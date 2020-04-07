#pragma once
#include "cfloat.h"
#include "typedefs.h"


#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/ref.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include <iostream>


constexpr double pi = 3.14159265358979323846;
constexpr cfloat I = { 0, 1 };
constexpr double sqrt2 = 1.41421356237309504880;
constexpr double sqrt2inv = 0.707106781186547524401;
constexpr double mu_0 = 1.25663706e-6; //vacuumm permeability
constexpr double mu_B = 927.4009994e-26; // J/T, Bohr magneton
constexpr double mu_N = 5.050783699e-27; // J/T, nuclear magneton
constexpr double hbar = 1.054571800e-34; // Js
constexpr double hplanck = 6.626070040e-34; // Js
constexpr double e = 1.6021766288e-19; // C
constexpr double m_e = 9.10938356e-31; // kg
constexpr double Celsius = 273.15; // K
constexpr double k_B = 1.38064852e-23; // J/K
constexpr double gfree = 2.00231930436182;
constexpr double y14N = 1.933e7;

constexpr double sqrt3 = 1.73205080757;
constexpr double sqrt6 = 2.44948974278;
constexpr double sqrt8 = 2.82842712475;
constexpr double sqrt38 = 0.61237243569;