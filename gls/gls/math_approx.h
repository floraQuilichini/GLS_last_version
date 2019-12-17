#pragma once
#include <cmath>
#include <algorithm>
#include <iostream>
#include "Eigen/Eigen"
#include "base.h"

class tanh_lookup_table
{
private:
	Scalar min_abscisse_;
	Scalar max_abscisse_;
	int nb_samples_;
	Scalar step_;
	std::vector<Scalar> lookup_table_;
public :
	tanh_lookup_table(Scalar min_abscisse, Scalar max_abscisse, int nb_samples);
	void read_lookup_table(const Eigen::ArrayXd& abscisses, Eigen::ArrayXd& tanh_values);

};
