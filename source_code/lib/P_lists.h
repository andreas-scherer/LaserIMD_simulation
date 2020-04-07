#pragma once
#include "typedefs.h"
#include "constants.h"
#include "functions.h"
#include "spin_system.h"
#include <filesystem>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/ref.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>




template<class Spin_System>
struct P_list_one_rot {

	Spin_System spin_system;
	double t_step;
	std::vector<double> theta;
	std::vector<double> phi;
	std::vector<typename Spin_System::operator_type> P_list;


	P_list_one_rot(Spin_System spin_system_, double t_step_, std::vector<double> theta_, std::vector<double> phi_) :
		spin_system(spin_system_),
		t_step(t_step_),
		theta(theta_),
		phi(phi_)
	{
		

		for (int i = 0; i < phi.size(); ++i) {
			P_list.push_back(typename Spin_System::operator_type::Identity());
		}

	}




	void create_P_list() {


		typename Spin_System::operator_type H;
		typename Spin_System::operator_type P;



	
		for (int k = 0; k < phi.size(); ++k) {

			std::vector<double> angles = { 0,theta[k],phi[k],0,theta[k],phi[k], 0,theta[k],phi[k] };
			H = spin_system.get_H(angles);
			
			P = (-I * H*t_step).exp().eval();
			P_list[k] = P;
		}

	}




	void copy_P_list(P_list_one_rot P1) {

		t_step = P1.t_step;

		for (int i = 0; i < phi.size(); ++i) {
			P_list[i] = P1.P_list[i];
		}
			
		
	}



	typename Spin_System::operator_type get_P_list(int i) {
		return P_list[i];
	}

	void reset() {
		t_step = 0.0;
		for (int i = 0; i < phi.size(); ++i) {
			P_list[i] = typename Spin_System::operator_type::Identity();
		}
	}

	void multiply(P_list_one_rot P1, bool forward) {

		if (forward) {
			t_step = t_step + P1.t_step;

			for (int i = 0; i < phi.size(); ++i) {
				P_list[i] = P_list[i] * P1.P_list[i];
			}


		}
		else {
			t_step = t_step - P1.t_step;

			for (int i = 0; i < phi.size(); ++i) {
				P_list[i] = P_list[i] * P1.P_list[i].adjoint().eval();
			}
		}



	}


};






template<class Spin_System>
struct P_list_three_rot {

	Spin_System spin_system;
	double t_step;
	std::vector<double> theta_obs;
	std::vector<double> phi_obs;
	std::vector<double> theta_T;
	std::vector<double> phi_T;
	std::vector<double> theta_dipolar;
	std::vector<double> phi_dipolar;
	std::vector<std::vector<std::vector<typename Spin_System::operator_type>>> P_list;


	P_list_three_rot(Spin_System spin_system_, double t_step_, std::vector<double> theta_obs_, std::vector<double> phi_obs_, std::vector<double> theta_T_, std::vector<double> phi_T_, std::vector<double> theta_dipolar_, std::vector<double> phi_dipolar_) :
		spin_system(spin_system_),
		t_step(t_step_),
		theta_obs(theta_obs_),
		phi_obs(phi_obs_),
		theta_T(theta_T_),
		phi_T(phi_T_),
		theta_dipolar(theta_dipolar_),
		phi_dipolar(phi_dipolar_)
	{
		std::vector<typename Spin_System::operator_type> dipolar;
		std::vector<std::vector<typename Spin_System::operator_type>> T;

		for (int k = 0; k < phi_dipolar.size(); ++k) {
			dipolar.push_back(typename Spin_System::operator_type::Identity());
		}

		for (int m = 0; m < phi_T.size(); ++m) {
			T.push_back(dipolar);
		}

		for (int i = 0; i < phi_obs.size(); ++i) {
			P_list.push_back(T);
		}
		
	}




	void create_P_list() {

		
		typename Spin_System::operator_type H;
		typename Spin_System::operator_type P;


		
		for (int i = 0; i < phi_obs.size(); ++i) {
			for (int m = 0; m <  phi_T.size(); ++m) {
				for (int k = 0; k < phi_dipolar.size(); ++k) {

					std::vector<double> angles = { 0,theta_obs[i],phi_obs[i], 0,theta_T[m],phi_T[m], 0,theta_dipolar[k],phi_dipolar[k] };
					H = spin_system.get_H(angles);
					P = (-I * H*t_step).exp().eval();
					P_list[i][m][k]=P;
				}
				
			}

		}
	} 

	


	void copy_P_list(P_list_three_rot P1) {

		t_step = P1.t_step;

		for (int i = 0; i < phi_obs.size(); ++i) {
			for (int m = 0; m < phi_T.size(); ++m) {
				for (int k = 0; k < phi_dipolar.size(); ++k) {

					P_list[i][m][k] = P1.P_list[i][m][k];

				}
			}
		}
	}



	typename Spin_System::operator_type get_P_list(int i, int k, int m) {
		return P_list[i][k][m];
	}

	void reset() {
		t_step = 0.0;
		for (int i = 0; i < phi_obs.size(); ++i) {
			for (int m = 0; m < phi_T.size(); ++m) {
				for (int k = 0; k < phi_dipolar.size(); ++k) {

					P_list[i][m][k] = typename Spin_System::operator_type::Identity();

				}
			}
		}
	}



	void multiply(P_list_three_rot P1, bool forward) {

		if (forward) {
			t_step = t_step + P1.t_step;

			for (int i = 0; i < phi_obs.size(); ++i) {
				for (int m = 0; m < phi_T.size(); ++m) {
					for (int k = 0; k < phi_dipolar.size(); ++k) {

						P_list[i][m][k] = P_list[i][m][k]*P1.P_list[i][m][k];

					}

				}

			}


		}
		else {
			t_step = t_step - P1.t_step;

			for (int i = 0; i < phi_obs.size(); ++i) {
				for (int m = 0; m < phi_T.size(); ++m) {
					for (int k = 0; k < phi_dipolar.size(); ++k) {

						P_list[i][m][k] = P_list[i][m][k] * P1.P_list[i][m][k].adjoint().eval();

					}

				}

			}
		}

		

	}


};



