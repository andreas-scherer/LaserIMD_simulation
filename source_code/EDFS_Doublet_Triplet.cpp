#include "Lib/spin_system.h"
#include "Lib/log_trace.h"
#include "Lib/simulate.h"
#include "Lib/grid.h"
#include "Lib/constants.h"
#include "Lib/load_parameters.h"

#include <omp.h>

#include <iostream>
#include <filesystem>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Eigenvalues>
#include <vector>
#include <fstream>
#include <string>
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/ref.hpp>





namespace fs = std::experimental::filesystem;

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



	std::vector<double> phi;
	std::vector<double> theta;
	std::vector<double> weights;


	std::vector<double> B;
	double t_step;



	Spin_System_Doublet_Triplet spin_time;
	Spin_System_Doublet_Triplet spin_pulse;
	std::vector<double> t;
	auto temp = boost::make_tuple(boost::ref(spin_time), boost::ref(spin_pulse), boost::ref(t_step), boost::ref(B), boost::ref(t));
	temp = load_parameters_Doublet_Triplet(argv[1]);




	auto temp2 = boost::make_tuple(boost::ref(phi), boost::ref(theta), boost::ref(weights));
	temp2 = grid(argv[3], spin_time.symmetrie_doublet, spin_time.knots_doublet);




	double const t0 = t[0];
	double const t1 = t[1];
	double const t2 = t[2];
	double const t3 = t[3];
	double const t4 = t[4];
	double const t5 = t[5];



	Log_stream log_stream(argv[2], "integral", t_step);



	for (auto B_temp : B) {


		spin_time.change_B(B_temp);
		spin_pulse.change_B(B_temp);





		log_stream.save_x(B_temp);

		for (int i = 0; i < phi.size(); ++i) {


			std::vector<double> angles = { 0,theta[i],phi[i], 0,theta[i],phi[i], 0,theta[i],phi[i],weights[i] };


			Spin_System_Doublet_Triplet::operator_type rho_start = spin_time.get_rho_equ(angles[0], angles[1], angles[2]);
			Spin_System_Doublet_Triplet::operator_type rho;

			rho = simulate(spin_pulse, t0, t1, angles, rho_start, Log_trace_none{});
			rho = simulate(spin_time, t1, t2, angles, rho, Log_trace_none{});
			rho = simulate(spin_pulse, t2, t3, angles, rho, Log_trace_none{});
			rho = simulate(spin_time, t3, t4, angles, rho, Log_trace_none{});
			rho = simulate(spin_time, t4, t5, angles, rho, log_stream);
		}

		log_stream.next();

	}



	log_stream.save_data();



}


