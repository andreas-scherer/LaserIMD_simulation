#include "Lib/spin_system.h"
#include "lib/log_trace.h"
#include "Lib/simulate.h"
#include "lib/grid.h"
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



	Spin_System_Doublet spin_time;
	Spin_System_Doublet spin_pulse;
	std::vector<double> t;
	auto temp = boost::make_tuple(boost::ref(spin_time), boost::ref(spin_pulse), boost::ref(t_step), boost::ref(B), boost::ref(t));
	temp = load_parameters_doublet(argv[1]);




	auto temp2 = boost::make_tuple(boost::ref(phi), boost::ref(theta), boost::ref(weights));
	temp2 = grid(argv[3], spin_time.symmetrie, spin_time.knots);




	double const t0 = t[0];
	double const t1 = t[1];
	double const t2 = t[2];
	double const t3 = t[3];
	double const t4 = t[4];
	double const t5 = t[5];



	Log_stream log_stream(argv[2], "transient", t_step);

	Spin_System_Doublet::operator_type rho_eq = spin_time.get_rho_equ();
	Spin_System_Doublet::operator_type rho;

	for (int i=0; i < phi.size(); ++i) {

		std::vector<double> angles = { 0,theta[i],phi[i],weights[i]};
		rho = simulate(spin_pulse, t0, t1, angles, rho_eq, Log_trace_none{});
		rho = simulate(spin_time, t1, t2, angles, rho, Log_trace_none{});
		rho = simulate(spin_pulse, t2, t3, angles, rho, Log_trace_none{});
		rho = simulate(spin_time, t3, t4, angles, rho, Log_trace_none{});
		rho = simulate(spin_time, t4, t5, angles, rho, log_stream);
		log_stream.reset();
	}
	
	log_stream.save_data();
	

}