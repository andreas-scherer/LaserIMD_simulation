#pragma once
#include <unsupported/Eigen/MatrixFunctions>
#include "constants.h"
#include "log_trace.h"
#include "functions.h"

#include <iosfwd>

template <class Spin_System, class Log_Stream>
typename Spin_System::operator_type simulate(Spin_System const & spin_system, double t_0, double t_end, std::vector<double> angles, typename Spin_System::operator_type const & initial_state, Log_Stream& log_output) {
	
	typename Spin_System::operator_type const H_0 = spin_system.get_H(angles);
	
	auto rho = initial_state;

	if (log_output.acq_mode!="none") {
		
		auto delta_t = log_output.t_step;
		
		typename Spin_System::operator_type P = (-I *H_0*delta_t).exp().eval();

		double t = t_0;
		for (t_0; t <= t_end; t += delta_t)
		{


			rho = P*rho*(P.adjoint()).eval();
		
			log_output.store<Spin_System>(t, rho, spin_system.Sx + spin_system.Tx +I*spin_system.Sy + I*spin_system.Ty, angles.back());
			
		}
		//do last step to the end

		P = (-I  * H_0*(t_end-t_0)).exp().eval();
		rho = P * initial_state*(P.adjoint().eval());
		
	} 	

	else	{
		
		

		auto P = ((-I * H_0*(t_end - t_0) ).exp()).eval();
		rho = P * rho*(P.adjoint().eval());
	


	}
	
	
	return rho;
}







template <class Spin_System, class Log_Stream>
typename Spin_System::operator_type simulate_P(Spin_System const & spin_system, typename Spin_System::operator_type const & P, double t_0, double t_end, double w, typename Spin_System::operator_type const & initial_state, Log_Stream& log_output) {

	//careful log_output.t_step must fit to t_step of propagator P

	auto rho = initial_state;
	
	if (log_output.acq_mode != "none") {

		auto delta_t = log_output.t_step;
		
		double t = t_0;
		for (t_0; t < t_end; t += delta_t)
		{

			rho = P * rho*(P.adjoint()).eval();
			log_output.store<Spin_System>(t, rho, spin_system.Sx + spin_system.Tx + I * spin_system.Sy + I * spin_system.Ty, w);   //     
		
		}
		

	}

	else {

		rho = P * rho*(P.adjoint().eval());
		
	}


	return rho;
}

