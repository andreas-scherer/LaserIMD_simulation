#include "Lib/spin_system.h"
#include "lib/log_trace.h"
#include "Lib/simulate.h"
#include "lib/grid.h"
#include "lib/constants.h"
#include "lib/P_lists.h"
#include "Lib/load_parameters.h"

#include <omp.h>

#include <iostream>
#include <filesystem>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <vector>
#include <fstream>
#include <string>
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/ref.hpp>
#include <math.h>

#include <chrono>


namespace fs = std::experimental::filesystem;

template<class P_list>
P_list create_8_ns_P_list(P_list const & P_list_2ns) {

	/*
	this function creates propagators that capture an 8ns step from propagators that capture a 2ns step
	*/

	//2ns
	P_list temp(P_list_2ns); 
	//4ns
	temp.multiply(temp,true);
	//8ns
	temp.multiply(temp, true);



	return temp;
}


template<class P_list>
P_list repeat_P_list(P_list const & P_list_step,double t_end) {

	/*
	this function repeats the propagators in P_list_step upto a time step t_end is captured
	*/

	P_list temp(P_list_step);
	temp.reset();


	while (temp.t_step < t_end ) {
		temp.multiply(P_list_step, true);
	}

	return temp;
}


template<class P_list>
P_list reverse_P_list(P_list const & P_list_start, P_list const & P_list_step, double t_end) {

	/*
	this function calculates propagators that capture a time step of the length t_end by reverting the
	propgatros in P_list_start with the propagators in P_list_step
	*/


	P_list temp(P_list_start);


	while ( temp.t_step > t_end) {
		temp.multiply(P_list_step, false);
		
	}

	return temp;
}




int main(int argc, char *argv[]) {

	try {
		if (argc < 4)
		{
			throw std::runtime_error("Not enough arguments given");
		}
	}
	catch (std::runtime_error &er) {
		std::cout << er.what() << "\n";
	}



	double t_step;
	std::vector<double> taus;
	Spin_System_Doublet_Triplet spin_time_2;
	Spin_System_Doublet_Triplet spin_pulse_2;
	std::vector<double> t;


	Spin_System_Doublet spin_time;
	Spin_System_Doublet spin_pulse;
		
	auto temp = boost::make_tuple(boost::ref(spin_time), boost::ref(spin_pulse), boost::ref(t_step), boost::ref(taus), boost::ref(t));
	temp = load_parameters_doublet(argv[1]);


	auto temp2 = boost::make_tuple(boost::ref(spin_time_2), boost::ref(spin_pulse_2), boost::ref(t_step), boost::ref(taus), boost::ref(t));
	temp2 = load_parameters_Doublet_Triplet(argv[1]);


	double const t0 = t[0];
	double const t1 = t[1];
	double const t2 = t[2];
	double const t3 = t[3];
	double const t4 = t[4];
	double const t5 = t[5];
	double t_pi_2 = t[1];

	double t_tau = taus[1] - taus[0];

	Log_stream log_stream(argv[2], "integral", t_step);

	std::vector<double> phi_doublet;
	std::vector<double> theta_doublet;
	std::vector<double> weights_doublet;
	auto temp_doublet = boost::make_tuple(boost::ref(phi_doublet), boost::ref(theta_doublet), boost::ref(weights_doublet));
	temp_doublet = grid(argv[3],spin_time_2.symmetrie_doublet, spin_time_2.knots_doublet);


	std::vector<double> phi_triplet;
	std::vector<double> theta_triplet;
	std::vector<double> weights_triplet;
	auto temp_triplet = boost::make_tuple(boost::ref(phi_triplet), boost::ref(theta_triplet), boost::ref(weights_triplet));
	temp_triplet = grid(argv[3], spin_time_2.symmetrie_triplet, spin_time_2.knots_triplet);

	std::vector<double> phi_dipolar;
	std::vector<double> theta_dipolar;
	std::vector<double> weights_dipolar;
	auto temp_dipolar = boost::make_tuple(boost::ref(phi_dipolar), boost::ref(theta_dipolar), boost::ref(weights_dipolar));
	temp_dipolar = grid(argv[3], spin_time_2.symmetrie_dipolar, spin_time_2.knots_dipolar);

	//the propagators for the Doublet system (before Laser excitation) are simulated in advance, they are updated later on capture the shift of the laser flash


	P_list_one_rot<Spin_System_Doublet> P_list_time_2ns(spin_time, t_step, theta_doublet, phi_doublet); P_list_time_2ns.create_P_list();
	P_list_one_rot<Spin_System_Doublet> P_list_time_8ns(create_8_ns_P_list(P_list_time_2ns));
	P_list_one_rot<Spin_System_Doublet> P_list_time_fid1(repeat_P_list(P_list_time_8ns, t2 - t1)); 
	P_list_one_rot<Spin_System_Doublet> P_list_time_fid2(repeat_P_list(P_list_time_8ns, t4 - t3));


	P_list_one_rot<Spin_System_Doublet> P_list_pulse_2ns(spin_pulse, t_step, theta_doublet, phi_doublet); P_list_pulse_2ns.create_P_list();
	P_list_one_rot<Spin_System_Doublet> P_list_pulse_8ns(create_8_ns_P_list(P_list_pulse_2ns)); 
	P_list_one_rot<Spin_System_Doublet>  P_list_pulse_pi_2(repeat_P_list(P_list_pulse_2ns, t_pi_2)); 
	P_list_one_rot<Spin_System_Doublet>  P_list_pulse_pi(P_list_pulse_pi_2);  P_list_pulse_pi.multiply(P_list_pulse_pi, true);


	P_list_one_rot<Spin_System_Doublet> P_list_dummy(P_list_pulse_2ns);
	Spin_System_Doublet::operator_type rho_eq = spin_time.get_rho_equ();
		
	Spin_System_Doublet::operator_type rho_doublet = Spin_System_Doublet::operator_type::Zero();
		

	std::vector<std::vector<Spin_System_Doublet::operator_type>> rho_doublets;

		


	for (auto tau : taus) {

		
		std::vector<Spin_System_Doublet::operator_type > rho_doublet_angles;
		
		log_stream.save_x(tau);
		
		std::cout << tau*1e9 << std::endl;

		if (tau < t0) {  //laser before first pulse


		}
		else if (tau >= t0 && tau < t1) {  //laser during first pulse

			if ((tau - t0) < t_tau) {
				P_list_dummy.copy_P_list(repeat_P_list(P_list_pulse_2ns, tau - t0));
			}
			else {
				P_list_dummy.multiply(P_list_pulse_8ns, true);
			}



		}
		else if (tau >= t1 && tau < t2) { //laser during first fid

			if ((tau - t1) < t_tau) {
				P_list_dummy.copy_P_list(repeat_P_list(P_list_time_2ns, tau - t1));
			}
			else {
				P_list_dummy.multiply(P_list_time_8ns, true);
			}
		}
		else if (tau >= t2 && tau < t3) { //laser during refocussing pulse

			if ((tau - t2) < t_tau) {
				P_list_dummy.copy_P_list(repeat_P_list(P_list_pulse_2ns, tau - t2));
			}
			else {
				P_list_dummy.multiply(P_list_pulse_8ns, true);
			}



		}
		else if (tau >= t3 && tau < t4) { //laser during second fid(echo)

			if ((tau - t3) < t_tau) {
				P_list_dummy.copy_P_list(repeat_P_list(P_list_time_2ns, tau - t3));
			}
			else {
				P_list_dummy.multiply(P_list_time_8ns, true);
			}
		}
		else {}



		for (int i = 0; i < phi_doublet.size(); ++i) {

				
			double w_doublet = weights_doublet[i];
			Spin_System_Doublet::operator_type rho;


			if (tau < t0) {  //laser before first pulse

				rho = rho_eq;

			}
			else if (tau >= t0 && tau < t1) {  //laser during first pulse
				rho = simulate_P(spin_pulse, P_list_dummy.P_list[i], t0, tau, w_doublet, rho_eq, Log_trace_none{});

			}
			else if (tau >= t1 && tau < t2) { //laser during first fid
				rho = simulate_P(spin_pulse, P_list_pulse_pi_2.P_list[i], t0, t1, w_doublet, rho_eq, Log_trace_none{});
				rho = simulate_P(spin_time, P_list_dummy.P_list[i], t1, tau, w_doublet, rho, Log_trace_none{});

			}
			else if (tau >= t2 && tau < t3) { //laser during refocussing pulse
				rho = simulate_P(spin_pulse, P_list_pulse_pi_2.P_list[i], t0, t1, w_doublet, rho_eq, Log_trace_none{});
				rho = simulate_P(spin_time, P_list_time_fid1.P_list[i], t1, t2, w_doublet, rho, Log_trace_none{});
				rho = simulate_P(spin_pulse, P_list_dummy.P_list[i], t2, tau, w_doublet, rho, Log_trace_none{});

			}
			else if (tau >= t3 && tau < t4) { //laser during second fid(echo)
				rho = simulate_P(spin_pulse, P_list_pulse_pi_2.P_list[i], t0, t1, w_doublet, rho_eq, Log_trace_none{});
				rho = simulate_P(spin_time, P_list_time_fid1.P_list[i], t1, t2, w_doublet, rho, Log_trace_none{});
				rho = simulate_P(spin_pulse, P_list_pulse_pi.P_list[i], t2, t3, w_doublet, rho, Log_trace_none{});
				rho = simulate_P(spin_time, P_list_dummy.P_list[i], t3, tau, w_doublet, rho, Log_trace_none{});

			}
			else if (tau >= t4 && tau < t5) { //laser during echo-detection
				rho = simulate_P(spin_pulse, P_list_pulse_pi_2.P_list[i], t0, t1, w_doublet, rho_eq, Log_trace_none{});
				rho = simulate_P(spin_time, P_list_time_fid1.P_list[i], t1, t2, w_doublet, rho, Log_trace_none{});
				rho = simulate_P(spin_pulse, P_list_pulse_pi.P_list[i], t2, t3, w_doublet, rho, Log_trace_none{});
				rho = simulate_P(spin_time, P_list_time_fid2.P_list[i], t3, t4, w_doublet, rho, Log_trace_none{});
				rho = simulate_P(spin_time, P_list_time_2ns.P_list[i], t4, tau, w_doublet, rho, log_stream);


			}
			else if (tau >= t5) { //laser after echo
				rho = simulate_P(spin_pulse, P_list_pulse_pi_2.P_list[i], t0, t1, w_doublet, rho_eq, Log_trace_none{});
				rho = simulate_P(spin_time, P_list_time_fid1.P_list[i], t1, t2, w_doublet, rho, Log_trace_none{});
				rho = simulate_P(spin_pulse, P_list_pulse_pi.P_list[i], t2, t3, w_doublet, rho, Log_trace_none{});
				rho = simulate_P(spin_time, P_list_time_fid2.P_list[i], t3, t4, w_doublet, rho, Log_trace_none{});
				rho = simulate_P(spin_time, P_list_time_2ns.P_list[i], t4, t5, w_doublet, rho, log_stream);
			}
			else {
			}
			rho_doublet_angles.push_back(rho);
		}

		rho_doublets.push_back(rho_doublet_angles);
		log_stream.next();
	}

	log_stream.reset();

		
	Eigen::Matrix3cd rho_T;


			
		
	spin_time_2.create_rho_T_list(theta_triplet, phi_triplet);
	
	bool first_part=false;

		
	//the propagators for the Doublet_Triplet system (after Laser excitation) are simulated in advance, they are updated later on capture the shift of the laser flash


	Spin_System_Doublet_Triplet::operator_type rho = Spin_System_Doublet_Triplet::operator_type::Zero();

	P_list_three_rot<Spin_System_Doublet_Triplet> P_list_full_time_2ns(spin_time_2, t_step, theta_doublet, phi_doublet, theta_triplet, phi_triplet, theta_dipolar, phi_dipolar); P_list_full_time_2ns.create_P_list();
	P_list_three_rot<Spin_System_Doublet_Triplet> P_list_full_time_8ns(create_8_ns_P_list(P_list_full_time_2ns));
	P_list_three_rot<Spin_System_Doublet_Triplet> P_list_full_time_fid1(repeat_P_list(P_list_full_time_8ns, t2 - t1));
	P_list_three_rot<Spin_System_Doublet_Triplet> P_list_full_time_fid2(repeat_P_list(P_list_full_time_8ns, t4 - t3));

	P_list_three_rot<Spin_System_Doublet_Triplet> P_list_full_pulse_2ns(spin_pulse_2, t_step, theta_doublet, phi_doublet, theta_triplet, phi_triplet, theta_dipolar, phi_dipolar); P_list_full_pulse_2ns.create_P_list();
	P_list_three_rot<Spin_System_Doublet_Triplet> P_list_full_pulse_8ns(create_8_ns_P_list(P_list_full_pulse_2ns));
	P_list_three_rot<Spin_System_Doublet_Triplet> P_list_full_pulse_pi_2(repeat_P_list(P_list_full_pulse_2ns, t_pi_2));
	P_list_three_rot<Spin_System_Doublet_Triplet> P_list_full_pulse_pi(P_list_full_pulse_pi_2);  P_list_full_pulse_pi.multiply(P_list_full_pulse_pi, true);

	P_list_three_rot<Spin_System_Doublet_Triplet> P_list_full_dummy(P_list_full_pulse_2ns);


	for (int q = 0; q < taus.size();++q) {
		double tau = taus[q];
		log_stream.save_x(tau);

		//show current tau on screen
		std::cout << tau*1e9 << "\n";
				
		for (int i = 0; i < phi_doublet.size(); ++i) {
					
			rho_doublet = rho_doublets[q][i];

			for (int k = 0; k < phi_dipolar.size(); ++k) {
				

				for (int m = 0; m < phi_triplet.size(); ++m) {

					rho_T = spin_time_2.get_rho_T_list(m);

					
					double w= weights_doublet[i] * weights_dipolar[k] * weights_triplet[m];
						
					if (tau < t0) {  //laser before first pulse

						if (!first_part) {
							rho = Eigen::kroneckerProduct(rho_eq, rho_T);
							rho = simulate_P(spin_pulse_2, P_list_full_pulse_pi_2.P_list[i][m][k], t0, t1, w, rho, Log_trace_none{});
							rho = simulate_P(spin_time_2, P_list_full_time_fid1.P_list[i][m][k], t1, t2, w, rho, Log_trace_none{});
							rho = simulate_P(spin_pulse_2, P_list_full_pulse_pi.P_list[i][m][k], t2, t3, w, rho, Log_trace_none{});
							rho = simulate_P(spin_time_2, P_list_full_time_fid2.P_list[i][m][k], t3, t4, w, rho, Log_trace_none{});
							rho = simulate_P(spin_time_2, P_list_full_time_2ns.P_list[i][m][k], t4, t5, w, rho, log_stream);
									
						}
						else {
							log_stream.copy_last_point();
							goto escape_point;
						}
								
					}
					else if (tau >= t0 && tau < t1) {  //laser during first pulse
							
						rho = Eigen::kroneckerProduct(rho_doublet, rho_T);
						rho = simulate_P(spin_pulse_2, P_list_full_dummy.P_list[i][m][k], tau, t1, w, rho, Log_trace_none{  });
						rho = simulate_P(spin_time_2, P_list_full_time_fid1.P_list[i][m][k], t1, t2, w, rho, Log_trace_none{  });
						rho = simulate_P(spin_pulse_2, P_list_full_pulse_pi.P_list[i][m][k], t2, t3, w, rho, Log_trace_none{  });
						rho = simulate_P(spin_time_2, P_list_full_time_fid2.P_list[i][m][k], t3, t4, w, rho, Log_trace_none{  });
						rho = simulate_P(spin_time_2, P_list_full_time_2ns.P_list[i][m][k], t4, t5, w, rho, log_stream);
					}
					else if (tau >= t1 && tau < t2) { //laser during first fid
								
						rho = Eigen::kroneckerProduct(rho_doublet, rho_T);
						rho = simulate_P(spin_time_2, P_list_full_dummy.P_list[i][m][k], tau, t2, w, rho, Log_trace_none{  });
						rho = simulate_P(spin_pulse_2, P_list_full_pulse_pi.P_list[i][m][k], t2, t3, w, rho, Log_trace_none{  });
						rho = simulate_P(spin_time_2, P_list_full_time_fid2.P_list[i][m][k], t3, t4, w, rho, Log_trace_none{  });
						rho = simulate_P(spin_time_2, P_list_full_time_2ns.P_list[i][m][k], t4, t5, w, rho, log_stream);

					}
					else if (tau >= t2 && tau < t3) { //laser during refocussing pulse
								
						rho = Eigen::kroneckerProduct(rho_doublet, rho_T);
						rho = simulate_P(spin_pulse_2, P_list_full_dummy.P_list[i][m][k], tau, t3, w, rho, Log_trace_none{  });
						rho = simulate_P(spin_time_2, P_list_full_time_fid2.P_list[i][m][k], t3, t4, w, rho, Log_trace_none{  });
						rho = simulate_P(spin_time_2, P_list_full_time_2ns.P_list[i][m][k], t4, t5, w, rho, log_stream);

					}
					else if (tau >= t3 && tau < t4) { //laser during second fid(echo)

						rho = Eigen::kroneckerProduct(rho_doublet, rho_T);
						rho = simulate_P(spin_time_2, P_list_full_dummy.P_list[i][m][k], tau, t4, w, rho, Log_trace_none{  });
						rho = simulate_P(spin_time_2, P_list_full_time_2ns.P_list[i][m][k], t4, t5, w, rho, log_stream);
						
					}
					else if (tau >= t4 && tau < t5) { //laser during echo-detection
						
						rho = Eigen::kroneckerProduct(rho_doublet, rho_T);
						rho = simulate_P(spin_time_2, P_list_full_time_2ns.P_list[i][m][k], tau, t5, w, rho, log_stream);
					}
					else if (tau >= t5) { //laser after echo

					}
					else {	}

		

				}//end T angle loop
			}//end dipolar angle loop
				
		}//end obs loop

		first_part = true;
		escape_point:;
		log_stream.next();

		if (tau < t0) {  //laser before first pulse
			if ((t0 - tau) < t_tau) {
				P_list_full_dummy.copy_P_list(reverse_P_list(P_list_full_pulse_pi_2, P_list_full_pulse_2ns, t1  - (tau + t_tau)));
			}
			else {

			}


		}
		else if (tau >= t0 && tau < t1) {  //laser during first pulse


			if ((t1 - tau) < t_tau) {
				P_list_full_dummy.copy_P_list(reverse_P_list(P_list_full_time_fid1, P_list_full_time_2ns, t2  - (tau + t_tau)));
			}
			else {
				P_list_full_dummy.multiply(P_list_full_pulse_8ns, false);
			}
		}
		else if (tau >= t1 && tau < t2) { //laser during first fid

			if ((t2 - tau) < t_tau) {
				P_list_full_dummy.copy_P_list(reverse_P_list(P_list_full_pulse_pi, P_list_full_pulse_2ns, t3 - (tau + t_tau)));
			}
			else {
				P_list_full_dummy.multiply(P_list_full_time_8ns, false);
			}

		}
		else if (tau >= t2 && tau < t3) { //laser during refocussing pulse

			if ((t3 - tau) < t_tau) {
				P_list_full_dummy.copy_P_list(reverse_P_list(P_list_full_time_fid2, P_list_full_time_2ns, t4  - (tau + t_tau)));
			}
			else {
				P_list_full_dummy.multiply(P_list_full_pulse_8ns, false);
			}

		}
		else if (tau >= t3 && tau < t4) { //laser during second fid(echo)

			if ((t4 - tau) < t_tau) {

			}
			else {
				P_list_full_dummy.multiply(P_list_full_time_8ns, false);
			}
		}
		else if (tau >= t4 && tau < t5) { //laser during echo-detection
		
		}
		else if (tau >= t5) { //laser after echo

		}
		else {}
			

	} //end tau loop		
	log_stream.reset();



	log_stream.save_data();
		
}


