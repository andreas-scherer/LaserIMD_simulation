#pragma once
#include "typedefs.h"
#include "constants.h"
#include "spin_system.h"

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <regex>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <Eigen/Eigenvalues>




namespace fs = std::experimental::filesystem;
//boost::tuple<Spin_System_Doublet, Spin_System_Doublet> 

boost::tuple<double, std::vector<double>, std::vector<double>> load_exp_par(std::string);

boost::tuple<Spin_System_Doublet, Spin_System_Doublet, double, std::vector<double>, std::vector<double>> load_parameters_doublet(std::string file) {

	std::ifstream fid(file);

	std::stringstream strStream;
	strStream << fid.rdbuf(); //read the file
	std::string str = strStream.str(); //str holds the content of the file

	

	// Load B
	std::smatch matches_B;
	std::regex B_regexp("B:[ ]+([\\d\\.]*)");
	double B;
	std::regex_search(str, matches_B, B_regexp);

	try {
		if (matches_B.size() == 0) {
			throw std::invalid_argument("no magnetic field given in proper format");
		}
			
		B = std::stod(matches_B[1]);
	

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}

	

	// Load g_D
	
	std::smatch matches_g_d;
	std::regex g_d_regexp("g_d:[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)");  //[ ]*([\\d\\.]*)

	Eigen::Matrix3d g_d;
	std::regex_search(str, matches_g_d, g_d_regexp);

	try {
		if (matches_g_d.size() == 0) {
			throw std::invalid_argument("no g value given in proper format");
		}

		g_d << std::stod(matches_g_d[1]), std::stod(matches_g_d[2]), std::stod(matches_g_d[3]),
			std::stod(matches_g_d[4]), std::stod(matches_g_d[5]), std::stod(matches_g_d[6]),
			std::stod(matches_g_d[7]), std::stod(matches_g_d[8]), std::stod(matches_g_d[9]);

		
	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}

	
	//Load flip period
	
	std::smatch matches_flip_period;
	std::regex flip_period_regexp("Flip period:[ ]+([\\d\\.]*)");
	double flip_period;
	std::regex_search(str, matches_flip_period, flip_period_regexp);

	try {
		if (matches_flip_period.size() == 0) {
			throw std::invalid_argument("no flip period given in proper format");
		}

		flip_period = std::stod(matches_flip_period[1])*1e-9;

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}

	
	double omega_rabi = 2.*pi / (flip_period);
	double B_nut = hbar * omega_rabi / (gfree*mu_B);



	//Load omega

	std::smatch matches_omega;
	std::regex omega_regexp("omega:[ ]+([\\d\\.]*)");
	double omega;
	std::regex_search(str, matches_omega, omega_regexp);

	try {
		if (matches_omega.size() == 0) {
			throw std::invalid_argument("no nutation frequency given in proper format");
		}

		omega = std::stod(matches_omega[1])*2*pi*1e9;

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	//Load temperature

	std::smatch matches_temperature;
	std::regex temperature_regexp("Temperature:[ ]+([\\d\\.]*)");
	double temperature;
	std::regex_search(str, matches_temperature, temperature_regexp);

	try {
		if (matches_temperature.size() == 0) {
			throw std::invalid_argument("no temperature given in proper format");
		}

		temperature = std::stod(matches_temperature[1]);

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}

	Eigen::Vector3d B_0;
	Eigen::Vector3d B_1;

	B_0 <<
		0,
		0,
		B;


	B_1 <<
		0,
		B_nut,
		B;


	Spin_System_Doublet spin_time(B_0, g_d, omega, temperature);
	Spin_System_Doublet spin_pulse(B_1, g_d, omega, temperature);


	// Load symmetrie
	std::smatch matches_symmetrie;
	std::regex symmetrie_regexp("Symmetrie doublet:[ ]+(Dinfh|Ci)");
	std::string symmetrie;
	std::regex_search(str, matches_symmetrie, symmetrie_regexp);

	try {
		if (matches_symmetrie.size() == 0) {
			throw std::invalid_argument("no Symmetrie given in proper format");
		}

		symmetrie=matches_symmetrie[1];


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	// Load knots
	std::smatch matches_knots;
	std::regex knots_regexp("Knots doublet:[ ]+([\\d]*)");
	int knots;
	std::regex_search(str, matches_knots, knots_regexp);

	try {
		if (matches_knots.size() == 0) {
			throw std::invalid_argument("no Knots given in proper format");
		}

		knots = std::stod(matches_knots[1]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	spin_time.set_symmetrie(symmetrie, knots);
	spin_pulse.set_symmetrie(symmetrie, knots);

	double delta_t;
	std::vector<double> t;
	std::vector<double> x;
	auto temp = boost::make_tuple(boost::ref(delta_t), boost::ref(x), boost::ref(t));
	temp = load_exp_par(str);

	auto res = boost::make_tuple(spin_time, spin_pulse, delta_t, x, t);
	return res;

}





boost::tuple<Spin_System_Doublet_Nuc_1, Spin_System_Doublet_Nuc_1, double, std::vector<double>, std::vector<double>> load_parameters_Doublet_Nuc_1(std::string file) {

	std::ifstream fid(file);

	std::stringstream strStream;
	strStream << fid.rdbuf(); //read the file
	std::string str = strStream.str(); //str holds the content of the file



									   // Load B
	std::smatch matches_B;
	std::regex B_regexp("B:[ ]+([\\d\\.]*)");
	double B;
	std::regex_search(str, matches_B, B_regexp);

	try {
		if (matches_B.size() == 0) {
			throw std::invalid_argument("no magnetic field given in proper format");
		}

		B = std::stod(matches_B[1]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}



	// Load g_D

	std::smatch matches_g_d;
	std::regex g_d_regexp("g_d:[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)");  //[ ]*([\\d\\.]*)

	Eigen::Matrix3d g_d;
	std::regex_search(str, matches_g_d, g_d_regexp);

	try {
		if (matches_g_d.size() == 0) {
			throw std::invalid_argument("no g value given in proper format");
		}

		g_d << std::stod(matches_g_d[1]), std::stod(matches_g_d[2]), std::stod(matches_g_d[3]),
			std::stod(matches_g_d[4]), std::stod(matches_g_d[5]), std::stod(matches_g_d[6]),
			std::stod(matches_g_d[7]), std::stod(matches_g_d[8]), std::stod(matches_g_d[9]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	//load HFI
	std::smatch matches_HFI;
	std::regex HFI_regexp("HFI:[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)");  //[ ]*([\\d\\.]*)

	Eigen::Matrix3d HFI;
	std::regex_search(str, matches_HFI, HFI_regexp);

	try {
		if (matches_HFI.size() == 0) {
			throw std::invalid_argument("no HFI value given in proper format");
		}

		HFI << std::stod(matches_HFI[1]), std::stod(matches_HFI[2]), std::stod(matches_HFI[3]),
			std::stod(matches_HFI[4]), std::stod(matches_HFI[5]), std::stod(matches_HFI[6]),
			std::stod(matches_HFI[7]), std::stod(matches_HFI[8]), std::stod(matches_HFI[9]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}

	HFI = pi * 2e6 * HFI;
	//Load flip period

	std::smatch matches_flip_period;
	std::regex flip_period_regexp("Flip period:[ ]+([\\d\\.]*)");
	double flip_period;
	std::regex_search(str, matches_flip_period, flip_period_regexp);

	try {
		if (matches_flip_period.size() == 0) {
			throw std::invalid_argument("no flip period given in proper format");
		}

		flip_period = std::stod(matches_flip_period[1])*1e-9;

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	double omega_rabi = 2.*pi / (flip_period);
	double B_nut = hbar * omega_rabi / (gfree*mu_B);



	//Load omega

	std::smatch matches_omega;
	std::regex omega_regexp("omega:[ ]+([\\d\\.]*)");
	double omega;
	std::regex_search(str, matches_omega, omega_regexp);

	try {
		if (matches_omega.size() == 0) {
			throw std::invalid_argument("no nutation frequency given in proper format");
		}

		omega = std::stod(matches_omega[1]) * 2 * pi*1e9;

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	//Load temperature

	std::smatch matches_temperature;
	std::regex temperature_regexp("Temperature:[ ]+([\\d\\.]*)");
	double temperature;
	std::regex_search(str, matches_temperature, temperature_regexp);

	try {
		if (matches_temperature.size() == 0) {
			throw std::invalid_argument("no temperature given in proper format");
		}

		temperature = std::stod(matches_temperature[1]);

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}

	Eigen::Vector3d B_0;
	Eigen::Vector3d B_1;

	B_0 <<
		0,
		0,
		B;


	B_1 <<
		0,
		B_nut,
		B;


	Spin_System_Doublet_Nuc_1 spin_time(B_0, g_d, y14N, HFI, omega, temperature);
	Spin_System_Doublet_Nuc_1 spin_pulse(B_1, g_d, y14N, HFI, omega, temperature);


	// Load symmetrie
	std::smatch matches_symmetrie;
	std::regex symmetrie_regexp("Symmetrie:[ ]+(Dinfh|Ci)");
	std::string symmetrie;
	std::regex_search(str, matches_symmetrie, symmetrie_regexp);

	try {
		if (matches_symmetrie.size() == 0) {
			throw std::invalid_argument("no Symmetrie given in proper format");
		}

		symmetrie = matches_symmetrie[1];


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	// Load knots
	std::smatch matches_knots;
	std::regex knots_regexp("Knots:[ ]+([\\d]*)");
	int knots;
	std::regex_search(str, matches_knots, knots_regexp);

	try {
		if (matches_knots.size() == 0) {
			throw std::invalid_argument("no Knots given in proper format");
		}

		knots = std::stod(matches_knots[1]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	spin_time.set_symmetrie(symmetrie, knots);
	spin_pulse.set_symmetrie(symmetrie, knots);

	double delta_t;
	std::vector<double> t;
	std::vector<double> x;
	auto temp = boost::make_tuple(boost::ref(delta_t), boost::ref(x), boost::ref(t));
	temp = load_exp_par(str);

	auto res = boost::make_tuple(spin_time, spin_pulse, delta_t, x, t);
	return res;

}








boost::tuple<Spin_System_Triplet, Spin_System_Triplet, double, std::vector<double>, std::vector<double>> load_parameters_Triplet(std::string file) {

	std::ifstream fid(file);

	std::stringstream strStream;
	strStream << fid.rdbuf(); //read the file
	std::string str = strStream.str(); //str holds the content of the file



									   // Load B
	std::smatch matches_B;
	std::regex B_regexp("B:[ ]+([\\d\\.]*)");
	double B;
	std::regex_search(str, matches_B, B_regexp);

	try {
		if (matches_B.size() == 0) {
			throw std::invalid_argument("no magnetic field given in proper format");
		}

		B = std::stod(matches_B[1]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}



	// Load g_t

	std::smatch matches_g_t;
	std::regex g_t_regexp("g_t:[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)");  //[ ]*([\\d\\.]*)

	Eigen::Matrix3d g_t;
	std::regex_search(str, matches_g_t, g_t_regexp);

	try {
		if (matches_g_t.size() == 0) {
			throw std::invalid_argument("no g value given in proper format");
		}

		g_t << std::stod(matches_g_t[1]), std::stod(matches_g_t[2]), std::stod(matches_g_t[3]),
			std::stod(matches_g_t[4]), std::stod(matches_g_t[5]), std::stod(matches_g_t[6]),
			std::stod(matches_g_t[7]), std::stod(matches_g_t[8]), std::stod(matches_g_t[9]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	//load D
	std::smatch matches_D;
	std::regex D_regexp("ZFS, D:[ ]+([-\\d\\.]*)");

	double D;
	std::regex_search(str, matches_D, D_regexp);

	try {
		if (matches_D.size() == 0) {
			throw std::invalid_argument("no D value given in proper format");
		}

		D = std::stod(matches_D[1]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}

	D = pi * 2e6 * D;



	//load E
	std::smatch matches_E;
	std::regex E_regexp("ZFS, E:[ ]+([-\\d\\.]*)");

	double E;
	std::regex_search(str, matches_E, E_regexp);

	try {
		if (matches_E.size() == 0) {
			throw std::invalid_argument("no E value given in proper format");
		}

		E = std::stod(matches_E[1]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}

	E = pi * 2e6 * E;


	//load Population
	std::smatch matches_P;
	std::regex P_regexp("Population:[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)");

	std::vector<double> P;
	std::regex_search(str, matches_P, P_regexp);

	try {
		if (matches_P.size() == 0) {
			throw std::invalid_argument("no P value given in proper format");
		}

		P.push_back(std::stod(matches_P[1]));
		P.push_back(std::stod(matches_P[2]));
		P.push_back(std::stod(matches_P[3]));


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}




	//Load flip period

	std::smatch matches_flip_period;
	std::regex flip_period_regexp("Flip period:[ ]+([\\d\\.]*)");
	double flip_period;
	std::regex_search(str, matches_flip_period, flip_period_regexp);

	try {
		if (matches_flip_period.size() == 0) {
			throw std::invalid_argument("no flip period given in proper format");
		}

		flip_period = std::stod(matches_flip_period[1])*1e-9;

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	double omega_rabi = 2.*pi / (flip_period);
	double B_nut = hbar * omega_rabi / (gfree*mu_B)/1.4142;



	//Load omega

	std::smatch matches_omega;
	std::regex omega_regexp("omega:[ ]+([\\d\\.]*)");
	double omega;
	std::regex_search(str, matches_omega, omega_regexp);

	try {
		if (matches_omega.size() == 0) {
			throw std::invalid_argument("no nutation frequency given in proper format");
		}

		omega = std::stod(matches_omega[1]) * 2 * pi*1e9;

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}



	Eigen::Vector3d B_0;
	Eigen::Vector3d B_1;

	B_0 <<
		0,
		0,
		B;


	B_1 <<
		0,
		B_nut,
		B;


	Spin_System_Triplet spin_time(B_0, g_t, D, E, P, omega);
	Spin_System_Triplet spin_pulse(B_1, g_t, D, E, P, omega);


	// Load symmetrie
	std::smatch matches_symmetrie;
	std::regex symmetrie_regexp("Symmetrie:[ ]+(Dinfh|Ci)");
	std::string symmetrie;
	std::regex_search(str, matches_symmetrie, symmetrie_regexp);

	try {
		if (matches_symmetrie.size() == 0) {
			throw std::invalid_argument("no Symmetrie given in proper format");
		}

		symmetrie = matches_symmetrie[1];


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	// Load knots
	std::smatch matches_knots;
	std::regex knots_regexp("Knots:[ ]+([\\d]*)");
	int knots;
	std::regex_search(str, matches_knots, knots_regexp);

	try {
		if (matches_knots.size() == 0) {
			throw std::invalid_argument("no Knots given in proper format");
		}

		knots = std::stod(matches_knots[1]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	spin_time.set_symmetrie(symmetrie, knots);
	spin_pulse.set_symmetrie(symmetrie, knots);



	double delta_t;
	std::vector<double> t;
	std::vector<double> x;
	auto temp = boost::make_tuple(boost::ref(delta_t), boost::ref(x), boost::ref(t));
	temp = load_exp_par(str);

	auto res = boost::make_tuple(spin_time, spin_pulse, delta_t, x, t);
	return res;

}


boost::tuple<Spin_System_Doublet_Triplet, Spin_System_Doublet_Triplet, double, std::vector<double>, std::vector<double>> load_parameters_Doublet_Triplet(std::string file) {

	std::ifstream fid(file);

	std::stringstream strStream;
	strStream << fid.rdbuf(); //read the file
	std::string str = strStream.str(); //str holds the content of the file



									   // Load B
	std::smatch matches_B;
	std::regex B_regexp("B:[ ]+([\\d\\.]*)");
	double B;
	std::regex_search(str, matches_B, B_regexp);

	try {
		if (matches_B.size() == 0) {
			throw std::invalid_argument("no magnetic field given in proper format");
		}

		B = std::stod(matches_B[1]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}




	// Load g_d

	std::smatch matches_g_d;
	std::regex g_d_regexp("g_d:[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)");  //[ ]*([\\d\\.]*)

	Eigen::Matrix3d g_d;
	std::regex_search(str, matches_g_d, g_d_regexp);

	try {
		if (matches_g_d.size() == 0) {
			throw std::invalid_argument("no g value given in proper format");
		}

		g_d << std::stod(matches_g_d[1]), std::stod(matches_g_d[2]), std::stod(matches_g_d[3]),
			std::stod(matches_g_d[4]), std::stod(matches_g_d[5]), std::stod(matches_g_d[6]),
			std::stod(matches_g_d[7]), std::stod(matches_g_d[8]), std::stod(matches_g_d[9]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}



	// Load g_t

	std::smatch matches_g_t;
	std::regex g_t_regexp("g_t:[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)");  //[ ]*([\\d\\.]*)

	Eigen::Matrix3d g_t;
	std::regex_search(str, matches_g_t, g_t_regexp);

	try {
		if (matches_g_t.size() == 0) {
			throw std::invalid_argument("no g value given in proper format");
		}

		g_t << std::stod(matches_g_t[1]), std::stod(matches_g_t[2]), std::stod(matches_g_t[3]),
			std::stod(matches_g_t[4]), std::stod(matches_g_t[5]), std::stod(matches_g_t[6]),
			std::stod(matches_g_t[7]), std::stod(matches_g_t[8]), std::stod(matches_g_t[9]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	//load D
	std::smatch matches_D;
	std::regex D_regexp("ZFS, D:[ ]+([-\\d\\.]*)");

	double D;
	std::regex_search(str, matches_D, D_regexp);

	try {
		if (matches_D.size() == 0) {
			throw std::invalid_argument("no D value given in proper format");
		}

		D = std::stod(matches_D[1]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}

	D = pi * 2e6 * D;



	//load E
	std::smatch matches_E;
	std::regex E_regexp("ZFS, E:[ ]+([-\\d\\.]*)");

	double E;
	std::regex_search(str, matches_E, E_regexp);

	try {
		if (matches_E.size() == 0) {
			throw std::invalid_argument("no E value given in proper format");
		}

		E = std::stod(matches_E[1]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}

	E = pi * 2e6 * E;


	//load Population
	std::smatch matches_P;
	std::regex P_regexp("Population:[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)[ ]+([\\d\\.]*)");

	std::vector<double> P;
	std::regex_search(str, matches_P, P_regexp);

	try {
		if (matches_P.size() == 0) {
			throw std::invalid_argument("no P value given in proper format");
		}

		P.push_back(std::stod(matches_P[1]));
		P.push_back(std::stod(matches_P[2]));
		P.push_back(std::stod(matches_P[3]));


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}




	//Load flip period

	std::smatch matches_flip_period;
	std::regex flip_period_regexp("Flip period:[ ]+([\\d\\.]*)");
	double flip_period;
	std::regex_search(str, matches_flip_period, flip_period_regexp);

	try {
		if (matches_flip_period.size() == 0) {
			throw std::invalid_argument("no flip period given in proper format");
		}

		flip_period = std::stod(matches_flip_period[1])*1e-9;

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	double omega_rabi = 2.*pi / (flip_period);
	double B_nut = hbar * omega_rabi / (gfree*mu_B);



	//Load omega

	std::smatch matches_omega;
	std::regex omega_regexp("omega:[ ]+([\\d\\.]*)");
	double omega;
	std::regex_search(str, matches_omega, omega_regexp);

	try {
		if (matches_omega.size() == 0) {
			throw std::invalid_argument("no nutation frequency given in proper format");
		}

		omega = std::stod(matches_omega[1]) * 2 * pi*1e9;

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	//Load temperature

	std::smatch matches_temperature;
	std::regex temperature_regexp("Temperature:[ ]+([\\d\\.]*)");
	double temperature;
	std::regex_search(str, matches_temperature, temperature_regexp);

	try {
		if (matches_temperature.size() == 0) {
			throw std::invalid_argument("no temperature given in proper format");
		}

		temperature = std::stod(matches_temperature[1]);

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}



	//Load t_yield

	std::smatch matches_t_yield;
	std::regex t_yield_regexp("T yield:[ ]+([\\d\\.]*)");
	double t_yield;
	std::regex_search(str, matches_t_yield, t_yield_regexp);

	try {
		if (matches_t_yield.size() == 0) {
			throw std::invalid_argument("no T yield given in proper format");
		}

		t_yield = std::stod(matches_t_yield[1]);

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}



	//Load distance

	std::smatch matches_distance;
	std::regex distance_regexp("Distance:[ ]+([\\d\\.]*)");
	double distance;
	std::regex_search(str, matches_distance, distance_regexp);

	try {
		if (matches_distance.size() == 0) {
			throw std::invalid_argument("no distance given in proper format");
		}

		distance = std::stod(matches_distance[1]);

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}

	double DD = (mu_0 / 2)*std::pow(mu_B, 2) * g_d.trace() / 3 * g_t.trace() / 3 / std::pow((distance*1e-9), 3) / hplanck;

	Eigen::Vector3d B_0;
	Eigen::Vector3d B_1;

	B_0 <<
		0,
		0,
		B;


	B_1 <<
		0,
		B_nut,
		B;


	Spin_System_Doublet_Triplet spin_time(B_0, g_d, g_t, DD, D, E, P, omega, t_yield, temperature);
	Spin_System_Doublet_Triplet spin_pulse(B_1, g_d, g_t, DD, D, E, P, omega, t_yield, temperature);


	// Load symmetrie doublet
	std::smatch matches_symmetrie_doublet;
	std::regex symmetrie_doublet_regexp("Symmetrie doublet:[ ]+(Dinfh|Ci)");
	std::string symmetrie_doublet;
	std::regex_search(str, matches_symmetrie_doublet, symmetrie_doublet_regexp);

	try {
		if (matches_symmetrie_doublet.size() == 0) {
			throw std::invalid_argument("no Symmetrie for the doublet given in proper format");
		}

		symmetrie_doublet = matches_symmetrie_doublet[1];


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	// Load knots doublet
	std::smatch matches_knots_doublet;
	std::regex knots_doublet_regexp("Knots doublet:[ ]+([\\d]*)");
	int knots_doublet;
	std::regex_search(str, matches_knots_doublet, knots_doublet_regexp);

	try {
		if (matches_knots_doublet.size() == 0) {
			throw std::invalid_argument("no knots for the doublet given in proper format");
		}

		knots_doublet = std::stod(matches_knots_doublet[1]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}



	// Load symmetrie triplet
	std::smatch matches_symmetrie_triplet;
	std::regex symmetrie_triplet_regexp("Symmetrie triplet:[ ]+(Dinfh|Ci)");
	std::string symmetrie_triplet;
	std::regex_search(str, matches_symmetrie_triplet, symmetrie_triplet_regexp);

	try {
		if (matches_symmetrie_triplet.size() == 0) {
			throw std::invalid_argument("no Symmetrie for the triplet given in proper format");
		}

		symmetrie_triplet = matches_symmetrie_triplet[1];


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	// Load knots triplet
	std::smatch matches_knots_triplet;
	std::regex knots_triplet_regexp("Knots triplet:[ ]+([\\d]*)");
	int knots_triplet;
	std::regex_search(str, matches_knots_triplet, knots_triplet_regexp);

	try {
		if (matches_knots_triplet.size() == 0) {
			throw std::invalid_argument("no knots for the triplet given in proper format");
		}

		knots_triplet = std::stod(matches_knots_triplet[1]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}




	// Load symmetrie dipolar
	std::smatch matches_symmetrie_dipolar;
	std::regex symmetrie_dipolar_regexp("Symmetrie dipolar:[ ]+(Dinfh|Ci)");
	std::string symmetrie_dipolar;
	std::regex_search(str, matches_symmetrie_dipolar, symmetrie_dipolar_regexp);

	try {
		if (matches_symmetrie_dipolar.size() == 0) {
			throw std::invalid_argument("no Symmetrie for the dipolar given in proper format");
		}

		symmetrie_dipolar = matches_symmetrie_dipolar[1];


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	// Load knots dipolar
	std::smatch matches_knots_dipolar;
	std::regex knots_dipolar_regexp("Knots dipolar:[ ]+([\\d]*)");
	int knots_dipolar;
	std::regex_search(str, matches_knots_dipolar, knots_dipolar_regexp);

	try {
		if (matches_knots_dipolar.size() == 0) {
			throw std::invalid_argument("no knots for the dipolar given in proper format");
		}

		knots_dipolar = std::stod(matches_knots_dipolar[1]);


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	spin_time.set_symmetrie(symmetrie_doublet, knots_doublet, symmetrie_triplet, knots_triplet, symmetrie_dipolar, knots_dipolar);
	spin_pulse.set_symmetrie(symmetrie_doublet, knots_doublet, symmetrie_triplet, knots_triplet, symmetrie_dipolar, knots_dipolar);

	double delta_t;
	std::vector<double> t;
	std::vector<double> x;
	auto temp = boost::make_tuple(boost::ref(delta_t), boost::ref(x), boost::ref(t));
	temp = load_exp_par(str);

	auto res = boost::make_tuple(spin_time, spin_pulse, delta_t, x, t);
	return res;

}


boost::tuple<double, std::vector<double>, std::vector<double>> load_exp_par(std::string str) {

	//Load flip period

	std::smatch matches_flip_period;
	std::regex flip_period_regexp("Flip period:[ ]+([\\d\\.]*)");
	double flip_period;
	std::regex_search(str, matches_flip_period, flip_period_regexp);

	try {
		if (matches_flip_period.size() == 0) {
			throw std::invalid_argument("no nutation frequency given in proper format");
		}

		flip_period = std::stod(matches_flip_period[1])*1e-9;

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";

	}

	// Load delta_t
	std::smatch matches_delta_t;
	std::regex delta_t_regexp("Delta t:[ ]+([\\d\\.]*)");
	double delta_t;
	std::regex_search(str, matches_delta_t, delta_t_regexp);

	try {
		if (matches_delta_t.size() == 0) {
			throw std::invalid_argument("no Delta t given in proper format");
		}

		delta_t = std::stod(matches_delta_t[1])*1e-9;


	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}



	//Load FID time

	std::smatch matches_fid_time;
	std::regex fid_time_regexp("FID time:[ ]+([\\d\\.]*)");
	double fid_time;
	std::regex_search(str, matches_fid_time, fid_time_regexp);

	try {
		if (matches_fid_time.size() == 0) {
			throw std::invalid_argument("no FID time given in proper format");
		}

		fid_time = std::stod(matches_fid_time[1])*1e-9;

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}



	//Load Integration window

	std::smatch matches_window;
	std::regex window_regexp("(Integration|Transient) window:[ ]+([\\d\\.]*)");
	double window;
	std::regex_search(str, matches_window, window_regexp);

	try {
		if (matches_window.size() == 0) {
			throw std::invalid_argument("no Integration or Transient window given in proper format");
		}

		window = std::stod(matches_window[2])*1e-9;

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	//Load Echo max

	std::smatch matches_echo_max;
	std::regex echo_max_regexp("Echo maximum:[ ]+([\\d\\.]*)");
	double echo_max;
	std::regex_search(str, matches_echo_max, echo_max_regexp);

	try {
		if (matches_echo_max.size() == 0) {
			throw std::invalid_argument("no Echo maximum given in proper format");
		}

		echo_max = std::stod(matches_echo_max[1])*1e-9;

	}
	catch (std::invalid_argument &er) {
		std::cout << "Exception occured" << "\n";
		std::cout << er.what() << "\n";
	}


	std::vector<double> x;

	if (matches_window[1].str() == "Integration") {
		// Load x_start
		std::smatch matches_x_start;
		std::regex x_start_regexp("(B|tau) start:[ ]+([-\\d\\.]*)");
		double x_start;
		std::regex_search(str, matches_x_start, x_start_regexp);

		try {
			if (matches_x_start.size() == 0) {
				throw std::invalid_argument("no x start given in proper format");
			}



			if (matches_x_start[1].str() == "B")
			{
				x_start = std::stod(matches_x_start[2]);
			}
			else if (matches_x_start[1].str() == "tau")
			{
				x_start = std::stod(matches_x_start[2])*1e-9;
			}


		}
		catch (std::invalid_argument &er) {
			std::cout << "Exception occured" << "\n";
			std::cout << er.what() << "\n";
		}



		// Load x_end
		std::smatch matches_x_end;
		std::regex x_end_regexp("(B|tau) end:[ ]+([-\\d\\.]*)");
		double x_end;
		std::regex_search(str, matches_x_end, x_end_regexp);

		try {
			if (matches_x_end.size() == 0) {
				throw std::invalid_argument("no x end given in proper format");
			}

			if (matches_x_end[1].str() == "B")
			{
				x_end = std::stod(matches_x_end[2]);
			}
			else if (matches_x_end[1].str() == "tau")
			{
				x_end = std::stod(matches_x_end[2])*1e-9;
			}


		}
		catch (std::invalid_argument &er) {
			std::cout << "Exception occured" << "\n";
			std::cout << er.what() << "\n";
		}



		// Load x_step
		std::smatch matches_x_step;
		std::regex x_step_regexp("(B|tau) step:[ ]+([-\\d\\.]*)");
		double x_step;
		std::regex_search(str, matches_x_step, x_step_regexp);

		try {
			if (matches_x_step.size() == 0) {
				throw std::invalid_argument("no x step given in proper format");
			}


			if (matches_x_step[1].str() == "B")
			{
				x_step = std::stod(matches_x_step[2]);
			}
			else if (matches_x_step[1].str() == "tau")
			{
				x_step = std::stod(matches_x_step[2])*1e-9;
			}


		}
		catch (std::invalid_argument &er) {
			std::cout << "Exception occured" << "\n";
			std::cout << er.what() << "\n";
		}


		x = arange(x_start, x_end, x_step);
	}

	std::vector<double> t;
	t.push_back(0);
	t.push_back(t[0] + flip_period / 4);
	t.push_back(t[1] + fid_time);
	t.push_back(t[2] + flip_period / 2);
	t.push_back(echo_max - window / 2);
	t.push_back(t[4] + window);

	auto res = boost::make_tuple(delta_t, x, t);
	return res;


}


