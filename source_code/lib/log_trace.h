#pragma once
#include <iostream>
#include <vector>
#include <cstddef>
#include <fstream>
#include <filesystem>

namespace fs = std::experimental::filesystem;

struct Log_trace_none {

	std::string file;
	std::string acq_mode;
	std::vector<double> x;
	std::vector<cfloat> y;
	int counter;
	double t_step;

	Log_trace_none()  :
		file(),
		acq_mode("none"),
		x(),
		y(),
		counter(0),
		t_step()	{}

	template <class Spin_System>
	void store(double t, typename Spin_System::operator_type const & rho, typename Spin_System::operator_type const & detc, double weight) {}

	template <class Spin_System>
	void store_additional(double t, typename Spin_System::operator_type const & rho, typename Spin_System::operator_type const &T_p_proj, typename Spin_System::operator_type const &T_0_proj, typename Spin_System::operator_type const &T_m_proj, double weight) {}
	

};


struct Log_stream{

	std::string file_name;
	std::string acq_mode;
	std::vector<double> x;
	std::vector<cfloat> y;
	std::vector<double> T_p;
	std::vector<double> T_0;
	std::vector<double> T_m;
	int counter;
	double add_total;
	double t_step;


	Log_stream(std::string file_name_, std::string acq_mode_, double t_step_) :
		file_name(file_name_),
		acq_mode(acq_mode_),
		x(),
		y(),
		T_p(),
		T_0(),
		T_m(),
		counter(0),
		add_total(0),
		t_step(t_step_)
	{
		
	
	}

	template <class Spin_System>
	void store(double t, typename Spin_System::operator_type const & rho, typename Spin_System::operator_type const &detc, double weight) {

		
		if (counter >= x.size()) {
			x.push_back(0);
			y.push_back(0);
		}
		
		if (acq_mode == "transient") {
			x[counter] = t;
		}

		y[counter] += weight * (detc*rho).trace();


		if (acq_mode == "transient") {
			++counter;
		}
		

	}


	void copy_last_point() {


		if (counter >= x.size()) {
			x.push_back(0);
			y.push_back(0);
		}


		y[counter] = y[counter - 1];



	}

	

	void save_x(double x_) {


		if (counter >= x.size()) {
			x.push_back(0);
			y.push_back(0);
		}

		x[counter] = x_;
	}

	void reset() { counter = 0; }

	void next() {
		++counter;
	}
	
	


	void save_data() {

		fs::path file(file_name);
		std::ofstream data_ascii(file, 'w');
		

		for (int i = 0; i < x.size(); ++i) {

			data_ascii << x[i] << "\t\t" << y[i].real() << "\t\t" << y[i].imag() << "\n";

		}


	}
		
};

