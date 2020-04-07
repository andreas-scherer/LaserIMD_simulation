#pragma once
#include "typedefs.h"
#include "constants.h"


#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/ref.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include <iostream>



boost::tuple<Eigen::MatrixXcd, Eigen::MatrixXcd, Eigen::MatrixXcd> sop(double i) {

	int n = static_cast<int>(2 * i + 1);
	
	Eigen::MatrixXcd Sx = Eigen::MatrixXcd::Zero(n, n);
	Eigen::MatrixXcd Sy = Eigen::MatrixXcd::Zero(n, n);
	Eigen::MatrixXcd Sz = Eigen::MatrixXcd::Zero(n, n);

	for (double x = -i, dim1 = n - 1; x <= i; ++x, --dim1) {

		for (double y = -i, dim2 = n - 1; y <= i; ++y, --dim2) {

			int d1 = static_cast<int>(dim1);
			int d2 = static_cast<int>(dim2);

			if (y == x + 1) {
				Sx(d1, d2) = 0.5*std::sqrt(i * (i + 1) - x * (x + 1));
				Sy(d1, d2) = I * 0.5*std::sqrt(i * (i + 1) - x * (x + 1));
			}

			if (y == x - 1) {
				Sx(d1, d2) = 0.5*std::sqrt(i * (i + 1) - x * (x - 1));
				Sy(d1, d2) = -I * 0.5*std::sqrt(i * (i + 1) - x * (x - 1));
			}

			if (y == x) {
				Sz(d1, d2) = x;

			}


		}

	}

	auto res = boost::make_tuple(Sx, Sy, Sz);
	return res;

}


template<class vec1, class vec2>
std::vector<vec2> all_terms(Eigen::Matrix3d a, std::vector<vec1> vec_op1, std::vector<vec2> vec_op2) {   //vec_op means vector operator can be R3 vector or Sx,Sy,Sz operators


	auto b = -a.trace() /sqrt3;
	auto T = -1 / sqrt3 * (vec_op1[0] * vec_op2[0] + vec_op1[1] * vec_op2[1] + vec_op1[2] * vec_op2[2]);


	auto Q = b * T;

	auto b_p2 = 0.5*(a(0, 0) - a(1, 1) - I * (a(0, 1) + a(1, 0)));
	auto b_p1 = -0.5*(a(0, 2) + a(2, 0) - I * (a(1, 2) + a(2, 1)));
	auto b_0 = 1 / sqrt6*(2 * a(2, 2) - (a(0, 0) + a(1, 1)));
	auto b_m1 = 0.5*(a(0, 2) + a(2, 0) + I * (a(1, 2) + a(2, 1)));
	auto b_m2 = 0.5*(a(0, 0) - a(1, 1) + I * (a(0, 1) + a(1, 0)));

	auto T_p2 = 1.0/2.0*(vec_op1[0] * vec_op2[0] - vec_op1[1] * vec_op2[1] - I * (vec_op1[0] * vec_op2[1] + vec_op1[1] * vec_op2[0]));
	auto T_p1 = -1.0 / 2.0*(vec_op1[0] * vec_op2[2] + vec_op1[2] * vec_op2[0] - I * (vec_op1[1] * vec_op2[2] + vec_op1[2] * vec_op2[1]));
	auto T_0 = 1 / sqrt6*(-vec_op1[0] * vec_op2[0] - vec_op1[1] * vec_op2[1] + 2 * vec_op1[2] * vec_op2[2]);
	auto T_m1 = 1.0 / 2.0*(vec_op1[0] * vec_op2[2] - vec_op1[2] * vec_op2[0] + I * (vec_op1[1] * vec_op2[2] + vec_op1[2] * vec_op2[1]));
	auto T_m2 = 1.0 / 2.0*(vec_op1[0] * vec_op2[0] - vec_op1[1] * vec_op2[1] + I * (vec_op1[0] * vec_op2[1] + vec_op1[1] * vec_op2[0]));


	auto Qm2m2 = b_m2 * T_m2;
	auto Qm2m1 = b_m1 * T_m2;
	auto Qm20 = b_0 * T_m2;
	auto Qm2p1 = b_p1 * T_m2;
	auto Qm2p2 = b_p2 * T_m2;

	auto Qm1m2 = b_m2 * T_m1;
	auto Qm1m1 = b_m1 * T_m1;
	auto Qm10 = b_0 * T_m1;
	auto Qm1p1 = b_p1 * T_m1;
	auto Qm1p2 = b_p2 * T_m1;

	auto Q0m2 = b_m2 * T_0;
	auto Q0m1 = b_m1 * T_0;
	auto Q00 = b_0 * T_0;
	auto Q0p1 = b_p1 * T_0;
	auto Q0p2 = b_p2 * T_0;

	auto Qp1m2 = b_m2 * T_p1;
	auto Qp1m1 = b_m1 * T_p1;
	auto Qp10 = b_0 * T_p1;
	auto Qp1p1 = b_p1 * T_p1;
	auto Qp1p2 = b_p2 * T_p1;

	auto Qp2m2 = b_m2 * T_p2;
	auto Qp2m1 = b_m1 * T_p2;
	auto Qp20 = b_0 * T_p2;
	auto Qp2p1 = b_p1 * T_p2;
	auto Qp2p2 = b_p2 * T_p2;

	std::vector<vec2> res = { Q,
		Qm2m2,Qm2m1,Qm20,Qm2p1,Qm2p2,
		Qm1m2,Qm1m1,Qm10,Qm1p1,Qm1p2,
		Q0m2,Q0m1,Q00,Q0p1,Q0p2,
		Qp1m2,Qp1m1,Qp10,Qp1p1,Qp1p2,
		Qp2m2,Qp2m1,Qp20,Qp2p1,Qp2p2 };

	return res;
}





template<class vec1, class vec2>
std::vector<vec2> sec_app(Eigen::Matrix3d a, std::vector<vec1> vec_op1, std::vector<vec2> vec_op2) {   //vec_op means vector operator can be R3 vector or Sx,Sy,Sz operators


	auto b = -a.trace() / sqrt3;
	auto T = -1 / sqrt3 * (vec_op1[0] * vec_op2[0] + vec_op1[1] * vec_op2[1] + vec_op1[2] * vec_op2[2]);


	auto Q = b * T;

	auto b_p2 = 0.5*(a(0, 0) - a(1, 1) - I * (a(0, 1) + a(1, 0)));
	auto b_p1 = -0.5*(a(0, 2) + a(2,0) - I * (a(1, 2) + a(2, 1)));
	auto b_0 = 1 / sqrt6*(2 * a(2, 2) - (a(0, 0) + a(1, 1)));
	auto b_m1 = 0.5*(a(0, 2) + a(2, 0) + I * (a(1, 2) + a(2, 1)));
	auto b_m2 = 0.5*(a(0, 0) - a(1, 1) + I * (a(0, 1) + a(1, 0)));

	auto T_p2 = vec2::Zero();
	auto T_p1 = vec2::Zero();
	auto T_0 = 1 / sqrt6*(-vec_op1[0] * vec_op2[0] - vec_op1[1] * vec_op2[1] + 2 * vec_op1[2] * vec_op2[2]);
	auto T_m1 = vec2::Zero();
	auto T_m2 = vec2::Zero();


	auto Qm2m2 = b_m2 * T_m2;
	auto Qm2m1 = b_m1 * T_m2;
	auto Qm20 = b_0 * T_m2;
	auto Qm2p1 = b_p1 * T_m2;
	auto Qm2p2 = b_p2 * T_m2;

	auto Qm1m2 = b_m2 * T_m1;
	auto Qm1m1 = b_m1 * T_m1;
	auto Qm10 = b_0 * T_m1;
	auto Qm1p1 = b_p1 * T_m1;
	auto Qm1p2 = b_p2 * T_m1;

	auto Q0m2 = b_m2 * T_0;
	auto Q0m1 = b_m1 * T_0;
	auto Q00 = b_0 * T_0;
	auto Q0p1 = b_p1 * T_0;
	auto Q0p2 = b_p2 * T_0;

	auto Qp1m2 = b_m2 * T_p1;
	auto Qp1m1 = b_m1 * T_p1;
	auto Qp10 = b_0 * T_p1;
	auto Qp1p1 = b_p1 * T_p1;
	auto Qp1p2 = b_p2 * T_p1;

	auto Qp2m2 = b_m2 * T_p2;
	auto Qp2m1 = b_m1 * T_p2;
	auto Qp20 = b_0 * T_p2;
	auto Qp2p1 = b_p1 * T_p2;
	auto Qp2p2 = b_p2 * T_p2;

	std::vector<vec2> res = { Q,
		Qm2m2,Qm2m1,Qm20,Qm2p1,Qm2p2,
		Qm1m2,Qm1m1,Qm10,Qm1p1,Qm1p2,
		Q0m2,Q0m1,Q00,Q0p1,Q0p2,
		Qp1m2,Qp1m1,Qp10,Qp1p1,Qp1p2,
		Qp2m2,Qp2m1,Qp20,Qp2p1,Qp2p2 };

	return res;
}


template<class vec1, class vec2>
std::vector<vec2> pseudo_sec_app(Eigen::Matrix3d a, std::vector<vec1> vec_op1, std::vector<vec2> vec_op2) {   //vec_op means vector operator can be R3 vector or Sx,Sy,Sz operators


	auto b = -a.trace() / sqrt3;
	auto T = -1 / sqrt3 * (vec_op1[2] * vec_op2[2]);


	auto Q = b * T;

	auto b_p2 = 0.5*(a(0, 0) - a(1, 1) - I * (a(0, 1) + a(1, 0)));
	auto b_p1 = -0.5*(a(0, 2) + a(2, 0) - I * (a(1, 2) + a(2, 1)));
	auto b_0 = 1 / sqrt6*(2 * a(2, 2) - (a(0, 0) + a(1, 1)));
	auto b_m1 = 0.5*(a(0, 2) + a(2, 0) + I * (a(1, 2) + a(2, 1)));
	auto b_m2 = 0.5*(a(0, 0) - a(1, 1) + I * (a(0, 1) + a(1, 0)));


	auto T_p2 = vec2::Zero();
	auto T_p1 = -0.5*(vec_op1[2] * (vec_op2[0] - I * vec_op2[1]));
	auto T_0 = 1 / sqrt6*(2 * vec_op1[2] * vec_op2[2]);
	auto T_m1 = 0.5*(vec_op1[2]*(vec_op2[0]+I*vec_op2[1]));
	auto T_m2 = vec2::Zero();


	auto Qm2m2 = b_m2 * T_m2;
	auto Qm2m1 = b_m1 * T_m2;
	auto Qm20 = b_0 * T_m2;
	auto Qm2p1 = b_p1 * T_m2;
	auto Qm2p2 = b_p2 * T_m2;

	auto Qm1m2 = b_m2 * T_m1;
	auto Qm1m1 = b_m1 * T_m1;
	auto Qm10 = b_0 * T_m1;
	auto Qm1p1 = b_p1 * T_m1;
	auto Qm1p2 = b_p2 * T_m1;

	auto Q0m2 = b_m2 * T_0;
	auto Q0m1 = b_m1 * T_0;
	auto Q00 = b_0 * T_0;
	auto Q0p1 = b_p1 * T_0;
	auto Q0p2 = b_p2 * T_0;

	auto Qp1m2 = b_m2 * T_p1;
	auto Qp1m1 = b_m1 * T_p1;
	auto Qp10 = b_0 * T_p1;
	auto Qp1p1 = b_p1 * T_p1;
	auto Qp1p2 = b_p2 * T_p1;

	auto Qp2m2 = b_m2 * T_p2;
	auto Qp2m1 = b_m1 * T_p2;
	auto Qp20 = b_0 * T_p2;
	auto Qp2p1 = b_p1 * T_p2;
	auto Qp2p2 = b_p2 * T_p2;

	


	std::vector<vec2> res = { Q,
		Qm2m2,Qm2m1,Qm20,Qm2p1,Qm2p2,
		Qm1m2,Qm1m1,Qm10,Qm1p1,Qm1p2,
		Q0m2,Q0m1,Q00,Q0p1,Q0p2,
		Qp1m2,Qp1m1,Qp10,Qp1p1,Qp1p2,
		Qp2m2,Qp2m1,Qp20,Qp2p1,Qp2p2 };




	return res;
}



template<class vec1, class vec2>
std::vector<vec2> weak_coupling(double omega, std::vector<vec1> vec_op1, std::vector<vec2> vec_op2) {   //vec_op means vector operator can be R3 vector or Sx,Sy,Sz operators


	auto b = omega;
	auto T = vec_op1[2] * vec_op2[2];


	auto Q = b * T;

	auto b_p2 = 0;
	auto b_p1 = 0;
	auto b_0 = 0;
	auto b_m1 = 0;
	auto b_m2 = 0;

	auto T_p2 = vec2::Zero();
	auto T_p1 = vec2::Zero();
	auto T_0 =  vec2::Zero();
	auto T_m1 = vec2::Zero();
	auto T_m2 = vec2::Zero();


	auto Qm2m2 = b_m2 * T_m2;
	auto Qm2m1 = b_m1 * T_m2;
	auto Qm20 = b_0 * T_m2;
	auto Qm2p1 = b_p1 * T_m2;
	auto Qm2p2 = b_p2 * T_m2;

	auto Qm1m2 = b_m2 * T_m1;
	auto Qm1m1 = b_m1 * T_m1;
	auto Qm10 = b_0 * T_m1;
	auto Qm1p1 = b_p1 * T_m1;
	auto Qm1p2 = b_p2 * T_m1;

	auto Q0m2 = b_m2 * T_0;
	auto Q0m1 = b_m1 * T_0;
	auto Q00 = b_0 * T_0;
	auto Q0p1 = b_p1 * T_0;
	auto Q0p2 = b_p2 * T_0;

	auto Qp1m2 = b_m2 * T_p1;
	auto Qp1m1 = b_m1 * T_p1;
	auto Qp10 = b_0 * T_p1;
	auto Qp1p1 = b_p1 * T_p1;
	auto Qp1p2 = b_p2 * T_p1;

	auto Qp2m2 = b_m2 * T_p2;
	auto Qp2m1 = b_m1 * T_p2;
	auto Qp20 = b_0 * T_p2;
	auto Qp2p1 = b_p1 * T_p2;
	auto Qp2p2 = b_p2 * T_p2;

	std::vector<vec2> res = { Q,
		Qm2m2,Qm2m1,Qm20,Qm2p1,Qm2p2,
		Qm1m2,Qm1m1,Qm10,Qm1p1,Qm1p2,
		Q0m2,Q0m1,Q00,Q0p1,Q0p2,
		Qp1m2,Qp1m1,Qp10,Qp1p1,Qp1p2,
		Qp2m2,Qp2m1,Qp20,Qp2p1,Qp2p2 };

	return res;
}




Eigen::MatrixXcd wignerm(int n, double alpha, double beta, double gamma) {

	Eigen::MatrixXcd Jx, Jy, Jz;
	auto temp = boost::make_tuple(boost::ref(Jx), boost::ref(Jy), boost::ref(Jz));
	temp = sop(n);

	Eigen::MatrixXcd d = ((-I * alpha *Jz).exp() * (-I * beta*Jy).exp() * (-I * gamma*Jz).exp()).eval();
	return d;

}



Matrix5cd wignerm2(double alpha, double beta, double gamma) {

	

	Matrix5cd d = Matrix5cd::Zero();
	Matrix5cd R = d;

	Eigen::Matrix<cfloat, 5, 1> a;
	a << std::exp(-2.0 *alpha* I),
		std::exp(-alpha * I),
		1,
		std::exp(alpha*I),
		std::exp(2.0 *alpha* I);

	Eigen::Matrix<cfloat, 1, 5> b;
	b << std::exp(-2.0 *gamma * I),
		std::exp(-gamma * I),
		1,
		std::exp(gamma*I),
		std::exp(2.0  * gamma * I);

	d(0, 0) = std::pow(1.0 + cos(beta), 2.0) / 4.0;
	d(0, 1) = -(1.0 + cos(beta))*sin(beta) / 2.0;
	d(0, 2) = sqrt38 * std::pow(sin(beta), 2);
	d(0, 3) = -(1.0 - cos(beta))*sin(beta) / 2.0;
	d(0, 4) = std::pow(1.0 - cos(beta), 2.0) / 4;

	d(1, 0) = (1.0 + cos(beta))*sin(beta) / 2.0;
	d(1, 1) = (cos(beta) - 1.0) / 2.0 + std::pow(cos(beta), 2.0);
	d(1, 2) = -sqrt38 * sin(2.0 * beta);
	d(1, 3) = (cos(beta) + 1.0) / 2.0 - std::pow(cos(beta), 2.0);
	d(1, 4) = -(1.0 - cos(beta))*sin(beta) / 2.0;

	d(2, 0) = sqrt38 * std::pow(sin(beta), 2.0);
	d(2, 1) = sqrt38 * sin(2 * beta);
	d(2, 2) = 0.5*(3.0 * std::pow(cos(beta), 2.0) - 1.0);
	d(2, 3) = -sqrt38 * sin(2.0 * beta);
	d(2, 4) = sqrt38 * std::pow(sin(beta), 2.0);

	d(3, 0) = (1.0 - cos(beta))*sin(beta) / 2.0;
	d(3, 1) = (cos(beta) + 1.0) / 2.0 - std::pow(cos(beta), 2.0);
	d(3, 2) = sqrt38 * sin(2 * beta);
	d(3, 3) = (cos(beta) - 1.0) / 2.0 + std::pow(cos(beta), 2.0);
	d(3, 4) = -(1.0 + cos(beta)) / 2.0 * sin(beta);

	d(4, 0) = std::pow(1.0 - cos(beta), 2.0) / 4.0;
	d(4, 1) = (1.0 - cos(beta))*sin(beta) / 2.0;
	d(4, 2) = sqrt38 * std::pow(sin(beta), 2.0);
	d(4, 3) = (1.0 + cos(beta))*sin(beta) / 2.0;
	d(4, 4) = 0.25*std::pow(1.0 + cos(beta), 2.0);

	R = (d.array()*(a*b).array()).matrix();
	return R;
}


template<class Spin_System>
std::vector<typename Spin_System::operator_type> plus(std::vector<typename Spin_System::operator_type> const& a, std::vector<typename Spin_System::operator_type> const &b) {

	if (a.size() != b.size()) {
		std::cerr << "Vectors must have the same length" << std::endl;
		
	}
	

	std::vector <typename Spin_System::operator_type> result( a.size() );

	for (int i = 0; i < a.size(); ++i) {
		result[i] = a[i] + b[i];
	}

	return result;

	

}


std::vector<double> linspace(double a, double b, int n) {
	std::vector<double> array;
	double step = (b - a) / (n - 1);

	while (a <= b) {
		array.push_back(a);
		a += step;           // could recode to better handle rounding errors
	}
	return array;
}


std::vector<double> arange(double a, double b, double step) {
	std::vector<double> array;

	while (a <= b) {
		array.push_back(a);
		a += step;           // could recode to better handle rounding errors
	}
	return array;
}


Eigen::MatrixXcd to_liouville(Eigen::MatrixXcd mat,int n)
{

	Eigen::MatrixXcd vec= Eigen::MatrixXcd::Zero(n*n,1);

	int k = 0;
	for (int j = 0; j < n; ++j)
	{
		for (int l = 0; l < n; ++l)
		{

			vec(k, 0) = mat(l, j);
			++k;
		}
	}
	return vec;
}


Eigen::MatrixXcd to_hilbert(Eigen::MatrixXcd vec, int n)
{
	Eigen::MatrixXcd mat = Eigen::MatrixXcd::Zero(n, n);

	int k = 0;
	for (int j = 0; j < n; ++j)
	{
		for (int l = 0; l < n; ++l)
		{

			mat(l, j)=vec(k, 0);
			++k;
		}
	}
	return mat;
}

