#pragma once
#include <cmath>
#include <algorithm>
#include <iostream>
#include "Eigen/Eigen"
#include "base.h"
#include "relative_scale.h"

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

Eigen::Matrix3d compute_pca_bbox_rescaling(PointMap& pointMap);
std::vector<VectorType> compute_regular_bbox_dims(Eigen::Matrix3Xd& data);