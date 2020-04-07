#pragma once
#include "typedefs.h"
#include "constants.h"

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <regex>
#include <filesystem>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/ref.hpp>
#include <boost/algorithm/string.hpp>



namespace fs = std::experimental::filesystem;

boost::tuple<std::vector<double>, std::vector<double>, std::vector<double>> grid(std::string folder, std::string symmetrie, int knots) {



	fs::path file (symmetrie+"_" +std::to_string(knots) + ".txt");
	fs::path path(folder);
	std::ifstream fid(path/file);
	std::vector<double> phi;
	std::vector<double> theta;
	std::vector<double> weights;
	if (!fid) {
		std::cerr << "File " << file << " cannot be opened\n";
	}
	else {
		std::string temp;
		std::getline(fid, temp);
		std::getline(fid, temp);
		do {
			
			
			std::vector<std::string> strs;
			boost::split(strs, temp, boost::is_any_of(";"));


			phi.push_back(std::stod(strs[0].c_str()));
			theta.push_back(std::stod(strs[1].c_str()));
			weights.push_back(std::stod(strs[2].c_str()));
			std::getline(fid, temp);
			

		} while (fid);

	}

	double sum = 0.0;
	for (double w : weights) {
		sum = sum + w;
	}

	for (double& w : weights) {
		w = w / sum;
	}
	
	auto result = boost::make_tuple(phi, theta, weights);
	return result;
}

