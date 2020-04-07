#pragma once
#include "typedefs.h"
#include "constants.h"
#include "functions.h"
#include <filesystem>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/ref.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>


struct Spin_System_Doublet {
	
	typedef Eigen::Matrix2cd operator_type;
	typedef Matrix4cd liouville_op_type;
	typedef Vector4cd lioville_state_type;

	Eigen::Vector3d B;
	Eigen::Matrix3d g_d;
	double omega_mw;
	double temperature;
	operator_type Sx;
	operator_type Sy;
	operator_type Sz;
	operator_type Tx;
	operator_type Ty;
	operator_type Tz;
	std::vector<operator_type> spher_tensor_op;
	std::string symmetrie;
	int knots;
	
	Spin_System_Doublet()
	{}


	Spin_System_Doublet(Eigen::Vector3d B_, Eigen::Matrix3d g_d_, double omega_mw_, double temperature_) :
		B(B_),
		g_d(g_d_),
		omega_mw(omega_mw_),
		temperature(temperature_)
	{
		auto temp = boost::make_tuple(boost::ref(Sx),boost::ref(Sy),  boost::ref(Sz));
		temp = sop(0.5);
		Tx = operator_type::Zero();
		Ty = operator_type::Zero();
		Tz = operator_type::Zero();


		
		std::vector<double> B_vec(3);
		B_vec[0] = B(0);
		B_vec[1] = B(1);
		B_vec[2] = B(2);
		spher_tensor_op = sec_app<double, operator_type>(mu_B/hbar* g_d, B_vec, std::vector<operator_type> {Sx, Sy, Sz});
		
	}


	operator_type get_H (std::vector<double> angles) const  {
		
		Matrix5cd d = wignerm2(angles[0], angles[1], angles[2]);
		
		
		operator_type H = spher_tensor_op[0] +
			d(0, 0)*spher_tensor_op[1] + d(0, 1)*spher_tensor_op[2] + d(0, 2)*spher_tensor_op[3] + d(0, 3)*spher_tensor_op[4] + d(0, 4)*spher_tensor_op[5] +
			d(1, 0)*spher_tensor_op[6] + d(1, 1)*spher_tensor_op[7] + d(1, 2)*spher_tensor_op[8] + d(1, 3)*spher_tensor_op[9] + d(1, 4)*spher_tensor_op[10] +
			d(2, 0)*spher_tensor_op[11] + d(2, 1)*spher_tensor_op[12] + d(2, 2)*spher_tensor_op[13] + d(2, 3)*spher_tensor_op[14] + d(2, 4)*spher_tensor_op[15] +
			d(3, 0)*spher_tensor_op[16] + d(3, 1)*spher_tensor_op[17] + d(3, 2)*spher_tensor_op[18] + d(3, 3)*spher_tensor_op[19] + d(3, 4)*spher_tensor_op[20] +
			d(4, 0)*spher_tensor_op[21] + d(4, 1)*spher_tensor_op[22] + d(4, 2)*spher_tensor_op[23] + d(4, 3)*spher_tensor_op[24] + d(4, 4)*spher_tensor_op[25]  -
			omega_mw * Sz;


		return H;
	
	
	}

	void set_symmetrie(std::string symmetrie_, int knots_)
	{
		symmetrie = symmetrie_;
		knots = knots_;
	}

	void change_B(double B_)
	{
		std::vector<double> B_vec(3);
		B_vec[0] = B(0);
		B_vec[1] = B(1);
		B_vec[2] = B_;
		spher_tensor_op = sec_app<double, operator_type>(mu_B / hbar * g_d, B_vec, std::vector<operator_type> {Sx, Sy, Sz});
	}

	operator_type R (operator_type rho) const  {
		return operator_type::Zero();
	}


	operator_type get_rho_equ()
	{
		Spin_System_Doublet::operator_type rho_eq = ((-hbar / k_B / temperature)*(omega_mw*Sz)).exp().eval(); rho_eq = rho_eq / rho_eq.trace();
		return rho_eq;
	}

	


};




struct Spin_System_Doublet_Nuc_1 {

	typedef Matrix6cd operator_type;

	Eigen::Vector3d B;
	Eigen::Matrix3d g_d;
	double yN;
	Eigen::Matrix3d HFI;
	double omega_mw;
	double temperature;
	operator_type Sx;
	operator_type Sy;
	operator_type Sz;
	operator_type Tx;
	operator_type Ty;
	operator_type Tz;
	operator_type Sx_nuc;
	operator_type Sy_nuc;
	operator_type Sz_nuc;
	std::vector<operator_type> spher_tensor_op;
	std::string symmetrie;
	int knots;

	Spin_System_Doublet_Nuc_1() {}

	Spin_System_Doublet_Nuc_1(Eigen::Vector3d B_, Eigen::Matrix3d g_d_, double yN_, Eigen::Matrix3d HFI_, double omega_mw_, double temperature_) :
		B(B_),
		g_d(g_d_),
		yN(yN_),
		HFI(HFI_),
		omega_mw(omega_mw_),
		temperature(temperature_)
	{
		
		
		//electron operators:
		Eigen::Matrix2cd Sx_temp;
		Eigen::Matrix2cd Sy_temp;
		Eigen::Matrix2cd Sz_temp;
		Eigen::Matrix3cd Unit2 = Eigen::Matrix3cd::Identity();
		auto temp = boost::make_tuple(boost::ref(Sx_temp), boost::ref(Sy_temp), boost::ref(Sz_temp));
		temp = sop(0.5);
		Sx = Eigen::kroneckerProduct(Sx_temp, Unit2);
		Sy = Eigen::kroneckerProduct(Sy_temp, Unit2);
		Sz = Eigen::kroneckerProduct(Sz_temp, Unit2);

		//nuclei operators
		Eigen::Matrix3cd Sx_temp2;
		Eigen::Matrix3cd Sy_temp2;
		Eigen::Matrix3cd Sz_temp2;
		Eigen::Matrix2cd Unit1 = Eigen::Matrix2cd::Identity();
		auto temp2 = boost::make_tuple(boost::ref(Sx_temp2), boost::ref(Sy_temp2), boost::ref(Sz_temp2));
		temp2 = sop(1);
		Sx_nuc = Eigen::kroneckerProduct(Unit1, Sx_temp2);
		Sy_nuc = Eigen::kroneckerProduct(Unit1, Sy_temp2);
		Sz_nuc = Eigen::kroneckerProduct(Unit1, Sz_temp2);

		
		Tx = operator_type::Zero();
		Ty = operator_type::Zero();
		Tz = operator_type::Zero();
		

		std::vector<double> B_vec(3);
		B_vec[0] = B(0);
		B_vec[1] = B(1);
		B_vec[2] = B(2);
		spher_tensor_op = sec_app<double, operator_type>(mu_B / hbar * g_d, B_vec, std::vector<operator_type> {Sx, Sy, Sz});
		spher_tensor_op = plus<Spin_System_Doublet_Nuc_1>(spher_tensor_op, pseudo_sec_app(HFI, std::vector<operator_type> {Sx, Sy, Sz}, std::vector<operator_type> {Sx_nuc, Sy_nuc, Sz_nuc}));
		spher_tensor_op[0] = spher_tensor_op[0] + yN * (B(0)*Sx_nuc + B(1)*Sy_nuc + B(2)*Sz_nuc);
	}


	void change_B(double B_)
	{
		std::vector<double> B_vec(3);
		B_vec[0] = B(0);
		B_vec[1] = B(1);
		B_vec[2] = B_;
		spher_tensor_op = sec_app<double, operator_type>(mu_B / hbar * g_d, B_vec, std::vector<operator_type> {Sx, Sy, Sz});
		spher_tensor_op = plus<Spin_System_Doublet_Nuc_1>(spher_tensor_op, pseudo_sec_app(HFI, std::vector<operator_type> {Sx, Sy, Sz}, std::vector<operator_type> {Sx_nuc, Sy_nuc, Sz_nuc}));
		spher_tensor_op[0] = spher_tensor_op[0] + yN * (B_vec[0]*Sx_nuc + B_vec[1]*Sy_nuc + B_vec[2]*Sz_nuc);
	}

	operator_type get_rho_equ()
	{
		operator_type rho_eq = ((-hbar / k_B / temperature)*(omega_mw*Sz)).exp().eval(); rho_eq = rho_eq / rho_eq.trace();
		return rho_eq;
	}


	operator_type get_H(std::vector<double> angles) const {

		Matrix5cd d = wignerm2(angles[0], angles[1], angles[2]);


		operator_type H = spher_tensor_op[0] +
			d(0, 0)*spher_tensor_op[1] + d(0, 1)*spher_tensor_op[2] + d(0, 2)*spher_tensor_op[3] + d(0, 3)*spher_tensor_op[4] + d(0, 4)*spher_tensor_op[5] +
			d(1, 0)*spher_tensor_op[6] + d(1, 1)*spher_tensor_op[7] + d(1, 2)*spher_tensor_op[8] + d(1, 3)*spher_tensor_op[9] + d(1, 4)*spher_tensor_op[10] +
			d(2, 0)*spher_tensor_op[11] + d(2, 1)*spher_tensor_op[12] + d(2, 2)*spher_tensor_op[13] + d(2, 3)*spher_tensor_op[14] + d(2, 4)*spher_tensor_op[15] +
			d(3, 0)*spher_tensor_op[16] + d(3, 1)*spher_tensor_op[17] + d(3, 2)*spher_tensor_op[18] + d(3, 3)*spher_tensor_op[19] + d(3, 4)*spher_tensor_op[20] +
			d(4, 0)*spher_tensor_op[21] + d(4, 1)*spher_tensor_op[22] + d(4, 2)*spher_tensor_op[23] + d(4, 3)*spher_tensor_op[24] + d(4, 4)*spher_tensor_op[25] -
			omega_mw * Sz;


		return H;


	}


	void set_symmetrie(std::string symmetrie_, int knots_)
	{
		symmetrie = symmetrie_;
		knots = knots_;
	}
	


};







struct Spin_System_Triplet {

	typedef Eigen::Matrix3cd operator_type;

	Eigen::Vector3d B;
	Eigen::Matrix3d g_t;
	double D;
	double E;
	std::vector<double> P;
	double omega_mw;
	operator_type Sx;
	operator_type Sy;
	operator_type Sz;
	operator_type Tx;
	operator_type Ty;
	operator_type Tz;
	std::vector<operator_type> spher_tensor_op;
	operator_type Bell;
	std::string symmetrie;
	int knots;

	Spin_System_Triplet() {}

	Spin_System_Triplet(Eigen::Vector3d B_, Eigen::Matrix3d g_t_, double D_, double E_, std::vector<double> P_, double omega_mw_) :
		B(B_),
		g_t(g_t_),
		D(D_),
		E(E_),
		P(P_),
		omega_mw(omega_mw_)
	{



		auto temp = boost::make_tuple(boost::ref(Tx), boost::ref(Ty), boost::ref(Tz));
		temp = sop(1);


		Sx = operator_type::Zero();
		Sy = operator_type::Zero();
		Sz = operator_type::Zero();

		Eigen::Matrix3d ZFS;
		ZFS << -1.0 / 3.0 * D + E, 0, 0,
			0, -1.0 / 3.0 * D - E, 0,
			0, 0, 2.0 / 3.0 * D;


		std::vector<double> B_vec(3);
		B_vec[0] = B(0);
		B_vec[1] = B(1);
		B_vec[2] = B(2);
		spher_tensor_op = sec_app<double, operator_type>(mu_B / hbar * g_t, B_vec, std::vector<operator_type> {Tx, Ty, Tz});
		spher_tensor_op = plus<Spin_System_Triplet>(spher_tensor_op, sec_app(ZFS, std::vector<operator_type> {Tx, Ty, Tz}, std::vector<operator_type> {Tx, Ty, Tz}));

		Bell << -1 / std::sqrt(2), I / std::sqrt(2), 0,
			0, 0, 1,
			1 / std::sqrt(2), I / std::sqrt(2), 0;

	
	}



	void change_B(double B_)
	{
		std::vector<double> B_vec(3);
		B_vec[0] = B(0);
		B_vec[1] = B(1);
		B_vec[2] = B_;

		Eigen::Matrix3d ZFS;
		ZFS << -1.0 / 3.0 * D + E, 0, 0,
			0, -1.0 / 3.0 * D - E, 0,
			0, 0, 2.0 / 3.0 * D;

		spher_tensor_op = sec_app<double, operator_type>(mu_B / hbar * g_t, B_vec, std::vector<operator_type> {Tx, Ty, Tz});
		spher_tensor_op = plus<Spin_System_Triplet>(spher_tensor_op, sec_app(ZFS, std::vector<operator_type> {Tx, Ty, Tz}, std::vector<operator_type> {Tx, Ty, Tz}));
	}

	



	operator_type get_H(std::vector<double> angles) const {

		Matrix5cd d = wignerm2(angles[0], angles[1], angles[2]);

		operator_type H = spher_tensor_op[0] +
			d(0, 0)*spher_tensor_op[1] + d(0, 1)*spher_tensor_op[2] + d(0, 2)*spher_tensor_op[3] + d(0, 3)*spher_tensor_op[4] + d(0, 4)*spher_tensor_op[5] +
			d(1, 0)*spher_tensor_op[6] + d(1, 1)*spher_tensor_op[7] + d(1, 2)*spher_tensor_op[8] + d(1, 3)*spher_tensor_op[9] + d(1, 4)*spher_tensor_op[10] +
			d(2, 0)*spher_tensor_op[11] + d(2, 1)*spher_tensor_op[12] + d(2, 2)*spher_tensor_op[13] + d(2, 3)*spher_tensor_op[14] + d(2, 4)*spher_tensor_op[15] +
			d(3, 0)*spher_tensor_op[16] + d(3, 1)*spher_tensor_op[17] + d(3, 2)*spher_tensor_op[18] + d(3, 3)*spher_tensor_op[19] + d(3, 4)*spher_tensor_op[20] +
			d(4, 0)*spher_tensor_op[21] + d(4, 1)*spher_tensor_op[22] + d(4, 2)*spher_tensor_op[23] + d(4, 3)*spher_tensor_op[24] + d(4, 4)*spher_tensor_op[25] -
			omega_mw * Tz;


		return H;


	}


	operator_type get_rho(double alpha, double beta, double gamma) const {

		

		operator_type rho;
		rho << P[0], 0, 0,
			0, P[1], 0,
			0, 0, P[2];
		
		operator_type R = wignerm(1, alpha, beta, gamma);
		
		rho = (R*Bell * rho*Bell.adjoint()*R.adjoint()).eval();

		rho << rho(0, 0), 0, 0,
			0, rho(1, 1), 0,
			0, 0, rho(2, 2);
		

		return rho;


	}

	void set_symmetrie(std::string symmetrie_, int knots_)
	{
		symmetrie = symmetrie_;
		knots = knots_;
	}


};







struct Spin_System_Doublet_Triplet {

	typedef Matrix6cd operator_type;
	typedef Matrix64cd liouville_op_type;
	typedef Vector64cd lioville_state_type;

	Eigen::Vector3d B;
	Eigen::Matrix3d g_d;
	Eigen::Matrix3d g_t;
	double DD;
	double D;
	double E;
	std::vector<double> P;
	double omega_mw;
	double t_yield;
	double temperature;
	operator_type Sx;
	operator_type Sy;
	operator_type Sz;
	operator_type Tx;
	operator_type Ty;
	operator_type Tz;
	Eigen::Matrix3cd Bell;
	std::vector<operator_type> spher_tensor_op_D;
	std::vector<operator_type> spher_tensor_op_T;
	std::vector<operator_type> spher_tensor_op_DD;
	std::vector<Eigen::Matrix3cd> rho_T_list;
	std::string symmetrie_doublet;
	int knots_doublet;
	std::string symmetrie_triplet;
	int knots_triplet;
	std::string symmetrie_dipolar;
	int knots_dipolar;

	Spin_System_Doublet_Triplet() {}

	Spin_System_Doublet_Triplet(Eigen::Vector3d B_, Eigen::Matrix3d g_d_, Eigen::Matrix3d g_t_, double DD_, double D_, double E_, std::vector<double> P_, double omega_mw_, double t_yield_, double temperature_) :
		B(B_),
		g_d(g_d_),
		g_t(g_t_),
		DD(DD_),
		D(D_),
		E(E_),
		P(P_),
		omega_mw(omega_mw_),
		t_yield(t_yield_),
		temperature(temperature_)
	{

		//doublet operators:
		Eigen::Matrix2cd Sx_temp;
		Eigen::Matrix2cd Sy_temp;
		Eigen::Matrix2cd Sz_temp;
		Eigen::Matrix3cd Unit2 = Eigen::Matrix3cd::Identity();
		auto temp = boost::make_tuple(boost::ref(Sx_temp), boost::ref(Sy_temp), boost::ref(Sz_temp));
		temp = sop(0.5);
		Sx = Eigen::kroneckerProduct(Sx_temp, Unit2);
		Sy = Eigen::kroneckerProduct(Sy_temp, Unit2);
		Sz = Eigen::kroneckerProduct(Sz_temp, Unit2);

		//triplet operators
		Eigen::Matrix3cd Tx_t;
		Eigen::Matrix3cd Ty_t;
		Eigen::Matrix3cd Tz_t;
		Eigen::Matrix2cd Unit1 = Eigen::Matrix2cd::Identity();


		
		auto temp2 = boost::make_tuple(boost::ref(Tx_t), boost::ref(Ty_t), boost::ref(Tz_t));
		temp2 = sop(1);



		Tx = Eigen::kroneckerProduct(Unit1, Tx_t);
		Ty = Eigen::kroneckerProduct(Unit1, Ty_t);
		Tz = Eigen::kroneckerProduct(Unit1, Tz_t);

		Eigen::Matrix3d DD_matrix;
		DD_matrix << -1, 0, 0,
			0, -1, 0,
			0, 0, 2;
		DD_matrix = DD * DD_matrix;

		Eigen::Matrix3d ZFS;
		ZFS << -1.0 / 3.0 * D + E, 0, 0,
			0, -1.0 / 3.0 * D - E, 0,
			0, 0, 2.0 / 3.0 * D;

		std::vector<double> B_vec(3);
		B_vec[0] = B(0);
		B_vec[1] = B(1);
		B_vec[2] = B(2);

		std::vector<operator_type> vec_op_S{ Sx, Sy, Sz };
		std::vector<operator_type> vec_op_T{ Tx, Ty, Tz };




		spher_tensor_op_D = sec_app<double, operator_type>(mu_B / hbar * g_d, B_vec, vec_op_S);
		spher_tensor_op_T = sec_app<double, operator_type>(mu_B / hbar * g_t, B_vec, vec_op_T);
		spher_tensor_op_T = plus<Spin_System_Doublet_Triplet>(spher_tensor_op_T, sec_app<operator_type, operator_type>(ZFS, vec_op_T, vec_op_T));
		spher_tensor_op_DD = sec_app<operator_type, operator_type>(DD_matrix, vec_op_S, vec_op_T);


		Bell << -1 / std::sqrt(2), I / std::sqrt(2), 0,
			0, 0, 1,
			1 / std::sqrt(2), I / std::sqrt(2), 0;
	}


	void change_distance(double DD_) {
		DD = DD_;
		Eigen::Matrix3d DD_matrix;
		DD_matrix << -1, 0, 0,
			0, -1, 0,
			0, 0, 2;
		DD_matrix = DD * DD_matrix;

		std::vector<operator_type> vec_op_S{ Sx, Sy, Sz };
		std::vector<operator_type> vec_op_T{ Tx, Ty, Tz };

		spher_tensor_op_DD = all_terms<operator_type, operator_type>(DD_matrix, vec_op_S, vec_op_T);
	}


	void change_B(double B_)
	{
		std::vector<double> B_vec(3);
		B_vec[0] = B(0);
		B_vec[1] = B(1);
		B_vec[2] = B_;

		Eigen::Matrix3d ZFS;
		ZFS << -1.0 / 3.0 * D + E, 0, 0,
			0, -1.0 / 3.0 * D - E, 0,
			0, 0, 2.0 / 3.0 * D;


		std::vector<operator_type> vec_op_S{ Sx, Sy, Sz };
		std::vector<operator_type> vec_op_T{ Tx, Ty, Tz };


		spher_tensor_op_D = sec_app<double, operator_type>(mu_B / hbar * g_d, B_vec, vec_op_S);
		spher_tensor_op_T = sec_app<double, operator_type>(mu_B / hbar * g_t, B_vec, vec_op_T);
		spher_tensor_op_T = plus<Spin_System_Doublet_Triplet>(spher_tensor_op_T, sec_app<operator_type, operator_type>(ZFS, vec_op_T, vec_op_T));

	}


	operator_type get_H(std::vector<double> angles) const {

		Matrix5cd D = wignerm2(angles[0], angles[1], angles[2]);
		Matrix5cd T = wignerm2(angles[3], angles[4], angles[5]);
		Matrix5cd DD = wignerm2(angles[6], angles[7], angles[8]);

		operator_type H = spher_tensor_op_D[0] +
			D(0, 0)*spher_tensor_op_D[1] + D(0, 1)*spher_tensor_op_D[2] + D(0, 2)*spher_tensor_op_D[3] + D(0, 3)*spher_tensor_op_D[4] + D(0, 4)*spher_tensor_op_D[5] +
			D(1, 0)*spher_tensor_op_D[6] + D(1, 1)*spher_tensor_op_D[7] + D(1, 2)*spher_tensor_op_D[8] + D(1, 3)*spher_tensor_op_D[9] + D(1, 4)*spher_tensor_op_D[10] +
			D(2, 0)*spher_tensor_op_D[11] + D(2, 1)*spher_tensor_op_D[12] + D(2, 2)*spher_tensor_op_D[13] + D(2, 3)*spher_tensor_op_D[14] + D(2, 4)*spher_tensor_op_D[15] +
			D(3, 0)*spher_tensor_op_D[16] + D(3, 1)*spher_tensor_op_D[17] + D(3, 2)*spher_tensor_op_D[18] + D(3, 3)*spher_tensor_op_D[19] + D(3, 4)*spher_tensor_op_D[20] +
			D(4, 0)*spher_tensor_op_D[21] + D(4, 1)*spher_tensor_op_D[22] + D(4, 2)*spher_tensor_op_D[23] + D(4, 3)*spher_tensor_op_D[24] + D(4, 4)*spher_tensor_op_D[25] -
			omega_mw * Sz;

		H = H + spher_tensor_op_T[0] +
			T(0, 0)*spher_tensor_op_T[1] + T(0, 1)*spher_tensor_op_T[2] + T(0, 2)*spher_tensor_op_T[3] + T(0, 3)*spher_tensor_op_T[4] + T(0, 4)*spher_tensor_op_T[5] +
			T(1, 0)*spher_tensor_op_T[6] + T(1, 1)*spher_tensor_op_T[7] + T(1, 2)*spher_tensor_op_T[8] + T(1, 3)*spher_tensor_op_T[9] + T(1, 4)*spher_tensor_op_T[10] +
			T(2, 0)*spher_tensor_op_T[11] + T(2, 1)*spher_tensor_op_T[12] + T(2, 2)*spher_tensor_op_T[13] + T(2, 3)*spher_tensor_op_T[14] + T(2, 4)*spher_tensor_op_T[15] +
			T(3, 0)*spher_tensor_op_T[16] + T(3, 1)*spher_tensor_op_T[17] + T(3, 2)*spher_tensor_op_T[18] + T(3, 3)*spher_tensor_op_T[19] + T(3, 4)*spher_tensor_op_T[20] +
			T(4, 0)*spher_tensor_op_T[21] + T(4, 1)*spher_tensor_op_T[22] + T(4, 2)*spher_tensor_op_T[23] + T(4, 3)*spher_tensor_op_T[24] + T(4, 4)*spher_tensor_op_T[25] -
			omega_mw * Tz;

		H = H + spher_tensor_op_DD[0] +
			DD(0, 0)*spher_tensor_op_DD[1] + DD(0, 1)*spher_tensor_op_DD[2] + DD(0, 2)*spher_tensor_op_DD[3] + DD(0, 3)*spher_tensor_op_DD[4] + DD(0, 4)*spher_tensor_op_DD[5] +
			DD(1, 0)*spher_tensor_op_DD[6] + DD(1, 1)*spher_tensor_op_DD[7] + DD(1, 2)*spher_tensor_op_DD[8] + DD(1, 3)*spher_tensor_op_DD[9] + DD(1, 4)*spher_tensor_op_DD[10] +
			DD(2, 0)*spher_tensor_op_DD[11] + DD(2, 1)*spher_tensor_op_DD[12] + DD(2, 2)*spher_tensor_op_DD[13] + DD(2, 3)*spher_tensor_op_DD[14] + DD(2, 4)*spher_tensor_op_DD[15] +
			DD(3, 0)*spher_tensor_op_DD[16] + DD(3, 1)*spher_tensor_op_DD[17] + DD(3, 2)*spher_tensor_op_DD[18] + DD(3, 3)*spher_tensor_op_DD[19] + DD(3, 4)*spher_tensor_op_DD[20] +
			DD(4, 0)*spher_tensor_op_DD[21] + DD(4, 1)*spher_tensor_op_DD[22] + DD(4, 2)*spher_tensor_op_DD[23] + DD(4, 3)*spher_tensor_op_DD[24] + DD(4, 4)*spher_tensor_op_DD[25];

		return H;

		
	}

	

	Eigen::Matrix3cd get_rho_T(double alpha, double beta, double gamma) const {


		Eigen::Matrix3cd rho_T;
		rho_T << P[0], 0, 0,
			0, P[1], 0,
			0, 0, P[2];
		Eigen::Matrix3cd R = wignerm(1, alpha, beta, gamma);

		rho_T = (R*Bell * rho_T*Bell.adjoint()*R.adjoint()).eval();
		Eigen::Matrix3cd rho = Eigen::Matrix3cd::Zero();
		rho(0, 0) = t_yield * rho_T(0, 0); rho(1, 1) = t_yield * rho_T(1, 1); rho(2, 2) = t_yield * rho_T(2, 2);

		return rho;
	}

	operator_type get_rho_full(Eigen::Matrix2cd rho_d, double alpha, double beta, double gamma) const {

		Eigen::Matrix3cd rho_T = get_rho_T(alpha, beta, gamma);
		operator_type rho = Eigen::kroneckerProduct(rho_d, rho_T);
		return rho;

	}


	operator_type get_rho_equ( double alpha, double beta, double gamma) const {

		Eigen::Matrix2cd Sx_temp;
		Eigen::Matrix2cd Sy_temp;
		Eigen::Matrix2cd Sz_temp;
		

		auto temp = boost::make_tuple(boost::ref(Sx_temp), boost::ref(Sy_temp), boost::ref(Sz_temp));
		temp = sop(0.5);
		Eigen::Matrix2cd rho_eq = ((-hbar / k_B / temperature)*(omega_mw*Sz_temp)).exp().eval(); rho_eq = rho_eq / rho_eq.trace();
		Eigen::Matrix3cd rho_T = get_rho_T(alpha, beta, gamma);
		operator_type rho = Eigen::kroneckerProduct(rho_eq, rho_T);
		return rho;

	}


	void set_symmetrie(std::string symmetrie_doublet_, int knots_doublet_, std::string symmetrie_triplet_, int knots_triplet_, std::string symmetrie_dipolar_, int knots_dipolar_)
	{
		symmetrie_doublet = symmetrie_doublet_;
		knots_doublet = knots_doublet_;
		symmetrie_triplet = symmetrie_triplet_;
		knots_triplet = knots_triplet_;
		symmetrie_dipolar = symmetrie_dipolar_;
		knots_dipolar = knots_dipolar_;
	}



	void create_rho_T_list(std::vector<double> theta_T, std::vector<double> phi_T) {
		
		Eigen::Matrix3cd rho_T;
		
		for (int m = 0; m < phi_T.size(); ++m) {

			rho_T = get_rho_T(0.0, theta_T[m], phi_T[m]);
			rho_T_list.push_back(rho_T);
		}

		
	} 

	Eigen::Matrix3cd get_rho_T_list(int m) {
		return rho_T_list[m];
	}



};





