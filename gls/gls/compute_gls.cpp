// gls.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


/*
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <cmath>
#include <algorithm>
#include <iostream>
#include "Patate/grenaille.h"
#include "Eigen/Eigen"
#include "relative_scale.h"
#include "base.h"
#include"IO.h"
#include "matching_prioritization.h"
#include "ransac_scheme.h"
#include <vector>
#include <tuple>
#include <fstream>
#include <string>
#include <chrono>

using namespace std;
using namespace Grenaille;


	/*---------------------- GLS computation --------------------------------*/

typedef Point::Scalar Scalar;
typedef Point::VectorType VectorType;
// Define related structure
typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightFunc;
typedef Basket<Point, WeightFunc, UnorientedSphereFit, GLSParam> Fit2;
typedef Basket<Point, WeightFunc, OrientedSphereFit, GLSParam> Fit1;
template<typename Fit>
std::tuple<Scalar, VectorType, Scalar, Scalar>  compute_gls_descriptor(Fit& _fit, vector<Point>& pointCloud, Scalar abs_scale, const Point& _p)
{

	// Set a weighting function instance
	_fit.setWeightFunc(WeightFunc(abs_scale));
	// Set the evaluation position
	_fit.init(_p.pos());  // choose point where you want to compute the descriptor

	// Iterate over samples and _fit the primitive
//#pragma omp parallel for -- not faster with parallel for
	for (int u = 0; u < pointCloud.size(); u++)
	{
		_fit.addNeighbor(pointCloud[u]);
	}
	//finalize fitting
	_fit.finalize();
	//Test if the fitting ended without errors
	if (_fit.isStable())
	{
		/*cout << "Center: [" << _fit.center().transpose() << "] ;  radius: " << _fit.radius() << endl;
		cout << "Pratt normalization"
			<< (_fit.applyPrattNorm() ? " is now done." : " has already been applied.") << endl;
		// Play with fitting output
		cout << "Fitted Sphere: " << endl
			<< "\t Tau : " << _fit.tau_normalized() << endl
			<< "\t Eta : " << _fit.eta_normalized().transpose() << endl
			<< "\t Kappa : " << _fit.kappa_normalized() << endl
			<< "\t Phi: " << _fit.fitness() << endl;*/

		std::tuple<Scalar, VectorType, Scalar, Scalar> gls_params = std::make_tuple(_fit.tau_normalized(), _fit.eta_normalized(), _fit.kappa_normalized(), _fit.fitness());
		return gls_params;

	}
	else
	{
		std::tuple<Scalar, VectorType, Scalar, Scalar> gls_params;
		return gls_params;
	}

}



typedef Point::Scalar Scalar;
typedef Point::VectorType VectorType;
// Define related structure
typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightFunc;
typedef Basket<Point, WeightFunc, OrientedSphereFit, GLSParam, OrientedSphereScaleDer, GLSDer, GLSGeomVar> DFit1;
template<typename DFit>
Scalar  compute_geometric_variation(DFit& _fit, vector<Point>& pointCloud, Scalar abs_scale, const Point& _p, VectorType W)
{
	Scalar geometric_variation = 0.0;

	// Set a weighting function instance
	_fit.setWeightFunc(WeightFunc(abs_scale));
	// Set the evaluation position
	_fit.init(_p.pos());  // choose point where you want to compute the descriptor

	// Iterate over samples and _fit the primitive
//#pragma omp parallel for
	for (int u = 0; u < pointCloud.size(); u++)
	{
		_fit.addNeighbor(pointCloud[u]);
	}
	//finalize fitting
	_fit.finalize();
	//Test if the fitting ended without errors

	if (_fit.isStable())
	{
		geometric_variation = _fit.geomVar(1.0, 1.0, 1.0);
		/*Scalar dtau = _fit.dtau_normalized().col(0)(0);
		Scalar deta = _fit.deta_normalized().col(0).norm();
		Scalar dkappa = _fit.dkappa_normalized().col(0)(0);
		Scalar geometric_variation = W[0]*dtau*dtau + W[1] * deta*deta + W[2] * dkappa*dkappa;*/
	}
	return geometric_variation;

	/*try {
		if (_fit.isStable())
		{
			geometric_variation = _fit.geomVar(1.0, 1.0, 1.0);
			return geometric_variation;
		}
		else
			throw std::string("fitting instable");
		
	}
	catch (std::string const& chaine)
	{
		std::cerr << chaine << std::endl;
	}
	*/
}


/*
int main(int argc, char** argv)
{

	if (argc < 5)
	{
		std::cout << "you must enter at least 4 arguments : the original point cloud input file, the down-sampled point cloud file where to compute the descriptors, the absolute scale of the point cloud, and the output file" << std::endl;
		std::cout << "If you want to compute GLS over a range of scales, you have to specify as arguments : the point cloud original input file, the down-sampled point cloud file where to compute the descriptors, the minimum scale, the maximum scale, the number of samples and the output file" << std::endl;
		return EXIT_FAILURE;
	}

	std::string pc_filename, pc_down_filename, output_filename;
	Scalar min_scale, max_scale, abs_scale, base;
	int nb_samples;
	bool flag_multiscale;
	bool flag_compute_geometric_variation = false;

	switch (argc)
	{
	case 8:
		flag_multiscale = true;
		pc_filename = argv[1];
		pc_down_filename = argv[2];
		min_scale = std::stod(argv[3]);
		max_scale = std::stod(argv[4]);
		nb_samples = std::stoi(argv[5]);
		output_filename = argv[6];
		flag_compute_geometric_variation = (std::stoi(argv[7]) > 0);
		break;
	case 7:
		flag_multiscale = true;
		pc_filename = argv[1];
		pc_down_filename = argv[2];
		min_scale = std::stod(argv[3]);
		max_scale = std::stod(argv[4]);
		nb_samples = std::stoi(argv[5]);
		output_filename = argv[6];
		break;
	case 6:
		flag_multiscale = false;
		pc_filename = argv[1];
		pc_down_filename = argv[2];
		min_scale = std::stod(argv[3]);
		output_filename = argv[4];
		max_scale = min_scale;
		nb_samples = 1;
		flag_compute_geometric_variation = (std::stoi(argv[7]) > 0);
		break;
	default:
		flag_multiscale = false;
		pc_filename = argv[1];
		pc_down_filename = argv[2];
		min_scale = std::stod(argv[3]);
		output_filename = argv[4];
		max_scale = min_scale;
		nb_samples = 1;
	}

	std::cout << "Input file:  " << pc_filename << std::endl;

	// get point cloud
	vector<Point> pointCloud, pointCloud_down;
	read_ply_file(pc_filename, pointCloud);
	read_ply_file(pc_down_filename, pointCloud_down);

	// set point where to compute descriptors
	//Point p = pointCloud[873];

	std::cout << "\n\n====================\nUnOrientedSphereFit:\n";
	//Fit2 fit2;
	Fit1 fit1;
	DFit1 dfit1;

	std::vector<std::pair<Point, std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>>>> points_gls_descriptors_over_scales;
	std::vector<std::pair<Point, std::vector<Scalar>>> points_geometric_variation_over_scales;

	if (flag_multiscale)
		base = pow(max_scale / min_scale, 1.0 / ((Scalar)nb_samples - 1.0));
	else
		base = 1.0;

//#pragma omp parallel for
	for (int j = 0; j < pointCloud_down.size(); j++)
	{
		Point p = pointCloud_down[j];
		std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>> point_gls_descriptors_over_scales;
		std::vector<Scalar> point_geometric_variation_over_scales;
//#pragma omp parallel for
		for (int i = 0; i < nb_samples; i++)
		{
			abs_scale = min_scale*pow(base, i);
			std::tuple<Scalar, VectorType, Scalar, Scalar> gls_params = compute_gls_descriptor(fit1, pointCloud, abs_scale, p);
			point_gls_descriptors_over_scales.push_back(gls_params);

			if (flag_compute_geometric_variation)
			{
				Scalar geom_var = compute_geometric_variation(dfit1, pointCloud, abs_scale, p, Eigen::Vector3d::Ones());
				point_geometric_variation_over_scales.push_back(geom_var);
				//std::cout << "geometric variation : " << geom_var << std::endl;
			}
		}
		points_gls_descriptors_over_scales.push_back(std::make_pair(p, point_gls_descriptors_over_scales));
		if (flag_compute_geometric_variation)
			points_geometric_variation_over_scales.push_back(std::make_pair(p, point_geometric_variation_over_scales));
	}

	// write output file
	int nb_points = pointCloud_down.size();
	if (flag_compute_geometric_variation)
		write_complete_profiles(nb_points, nb_samples, min_scale, max_scale, base, points_gls_descriptors_over_scales, output_filename, &points_geometric_variation_over_scales);
	else
		write_complete_profiles(nb_points, nb_samples, min_scale, max_scale, base, points_gls_descriptors_over_scales, output_filename);

	return 0;

}
*/

/*
int main(int argc, char** argv)
{

	if (argc < 4)
	{
		std::cout << "you must enter at least 3 arguments : the source descriptors file, the target descriptors file and the output file " << std::endl;
		return EXIT_FAILURE;
	}

	std::string descriptors_source_filename, descriptors_target_filename, output_filename;
	Scalar min_scale, max_scale, base;
	int nb_samples, nb_source_points, nb_target_points;
	bool read_geom_var = false;

	descriptors_source_filename = argv[1];
	descriptors_target_filename = argv[2];
	output_filename = argv[3];
	if(argc == 5)
		read_geom_var = (std::stoi(argv[4]) > 0);

	// get descriptors and headers 
	std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>> source_descriptors, target_descriptors;
	std::map<Point, std::vector<Scalar>, PointComp> source_geom_var, target_geom_var;
	if (read_geom_var)
	{
		read_descriptors_text_file(descriptors_source_filename, source_descriptors, nb_source_points, nb_samples, min_scale, max_scale, base, &source_geom_var);
		read_descriptors_text_file(descriptors_target_filename, target_descriptors, nb_target_points, nb_samples, min_scale, max_scale, base, &target_geom_var);
	}
	else
	{
		read_descriptors_text_file(descriptors_source_filename, source_descriptors, nb_source_points, nb_samples, min_scale, max_scale, base);
		read_descriptors_text_file(descriptors_target_filename, target_descriptors, nb_target_points, nb_samples, min_scale, max_scale, base);
	}

	// get relative scale estimation
	Scalar alpha = 1.0;
	VectorType w = Eigen::Vector3d::Ones();
	Scalar ratio = 0.1;

	std::pair<Scalar, std::vector<std::tuple<Point, Point, Scalar>>> relative_scale_and_matching_points = relative_scale_estimation(source_descriptors, target_descriptors, ratio, nb_samples, w, alpha);
	std::cout << "relative scale: " << relative_scale_and_matching_points.first << std::endl;

	if (read_geom_var)
	{
		pair_priority_queue queue;
		
		for (int i = 0; i < relative_scale_and_matching_points.second.size(); i++)
		{
			// get matching pair cost
			Point source = std::get<0>(relative_scale_and_matching_points.second[i]);
			Point target = std::get<1>(relative_scale_and_matching_points.second[i]);
			Scalar relative_scale = std::get<2>(relative_scale_and_matching_points.second[i]);

			std::map<Point, std::vector<Scalar>, PointComp>::iterator it_source = source_geom_var.find(source);
			std::map<Point, std::vector<Scalar>, PointComp>::iterator it_target = target_geom_var.find(target);

			if (it_source != source_geom_var.end() && it_target != target_geom_var.end())
			{
				// fill priority queue
				std::tuple<Point, Point, Scalar>  pair_cost = compute_pair_cost(it_source, it_target, nb_samples, alpha);
				queue.add_pair(std::make_tuple(source, target, relative_scale, std::get<2>(pair_cost)));
			}
			
		}
	
	}

	// write relative scale and matching points in file
	//std::string output_filename = "C:\\Registration\\test_multiscale\\output\\scales_matching_points.txt";
	write_relative_scale_and_matching_points(output_filename, relative_scale_and_matching_points);

	return 0;

}
*/


/*
int main(int argc, char** argv)
{

	if (argc < 4)
	{
		std::cout << "you must enter at least 3 arguments : the source descriptors file, the target descriptors file and the output file " << std::endl;
		return EXIT_FAILURE;
	}

	std::string descriptors_source_filename, descriptors_target_filename, output_filename;
	Scalar min_scale, max_scale, base;
	int nb_samples, nb_source_points, nb_target_points;

	descriptors_source_filename = argv[1];
	descriptors_target_filename = argv[2];
	output_filename = argv[3];

	// get descriptors and headers 
	std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>> source_descriptors, target_descriptors;
	read_descriptors_text_file(descriptors_source_filename, source_descriptors, nb_source_points, nb_samples, min_scale, max_scale, base);
	read_descriptors_text_file(descriptors_target_filename, target_descriptors, nb_target_points, nb_samples, min_scale, max_scale, base);

	// get matching pairs
	Scalar alpha = 1.0;
	VectorType w = Eigen::Vector3d::Ones();

	std::vector<std::pair<Point, Point>> matching_points = compute_matching_pairs_after_rescaling(source_descriptors, target_descriptors, w, alpha, nb_samples, true);

	// write matching points in file
	write_matching_points(output_filename, matching_points);

	return 0;

}
*/


int main(int argc, char** argv)
{

	if (argc < 4)
	{
		std::cout << "you must enter at least 3 arguments : the source descriptors file,  the target descriptors file, and the output file where to write the transform matrix" << std::endl;
		return EXIT_FAILURE;
	}

	std::string descriptors_source_filename, descriptors_target_filename, transform_filename;
	Scalar min_source_scale, max_source_scale, min_target_scale, max_target_scale, base;
	int nb_source_samples, nb_target_samples, nb_source_points, nb_target_points;

	descriptors_source_filename = argv[1];
	descriptors_target_filename = argv[2];
	transform_filename = argv[3];
	Scalar alpha = 4.0;

	// get descriptors and headers 
	std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>> source_descriptors_geom_var, target_descriptors_geom_var;
	read_descriptors_text_file(descriptors_source_filename, source_descriptors_geom_var, nb_source_points, nb_source_samples, min_source_scale, max_source_scale, base);
	read_descriptors_text_file(descriptors_target_filename, target_descriptors_geom_var, nb_target_points, nb_target_samples, min_target_scale, max_target_scale, base);
	
	// fill point maps
	PointMap map_index_profiles_cost_source, map_index_profiles_cost_target;
	map_index_profiles_cost_source.set_point_profiles_cost(source_descriptors_geom_var, nb_source_samples, alpha);
	map_index_profiles_cost_target.set_point_profiles_cost(target_descriptors_geom_var, nb_target_samples, alpha);

	// compute matching pairs with priority
	VectorType w = Eigen::Vector3d::Ones();
	Scalar ratio = 0.2;
	//std::vector<std::tuple<Point, Point, Scalar, Scalar>> three_closest_pairs = compute_3_closest_pairs(source_descriptors_geom_var, target_descriptors_geom_var, ratio, nb_source_samples, nb_target_samples, w, 4.0);
	std::vector<std::tuple<std::pair<int, int>, Scalar, Scalar, Scalar>> symmetric_pairs = compute_symmetric_pairs(map_index_profiles_cost_source, map_index_profiles_cost_target, ratio, 3, w, 4.0, false);
	
	// apply Ransac scheme
	int nb_iterations = 1000;
	Scalar max_err_scale = 0.2;
	Scalar max_err_reg = 0.8; 
	Scalar max_err_norm = 0.35;
	int nb_kpairs = 12;
	Scalar bbox_diag_ratio = 8.0;
	//pair_priority_queue queue(three_closest_pairs, min_source_scale, min_target_scale, base);
	pair_priority_queue queue(symmetric_pairs, &map_index_profiles_cost_source, &map_index_profiles_cost_target, min_source_scale, min_target_scale, base);
	//queue.display_queue();
	std::string output_filename = "C:\\Registration\\test_gls_algo\\matching_pairs\\source_matching_pairs.txt";
	std::string kpairs_filename = "C:\\Registration\\test_gls_algo\\matching_pairs\\kpairs_file.txt";
	//std::string tuple_filename = "C:\\Registration\\test_gls_algo\\matching_pairs\\tuple_file.txt";
	write_closest_matching_points(queue, &map_index_profiles_cost_source, &map_index_profiles_cost_target, output_filename, false);
	write_kpairs(kpairs_filename, queue, &map_index_profiles_cost_source, &map_index_profiles_cost_target, max_err_scale, nb_kpairs, bbox_diag_ratio);
	//write_tuples(tuple_filename, queue, &map_index_profiles_cost_source, &map_index_profiles_cost_target, max_err_scale);
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	RansacScheme ransac(queue, &map_index_profiles_cost_source, &map_index_profiles_cost_target);
	Eigen::Matrix4d transform = ransac.ransac_algorithm(nb_iterations, max_err_scale, max_err_reg, max_err_norm, bbox_diag_ratio);
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::chrono::nanoseconds ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0);
	std::cout << "It took " << ns.count() << " nanosecond(s) to run RANSAC" << std::endl;

	//std::cout << "transform : " << transform << std::endl;

	// write transform
	write_matrix_transform(transform, transform_filename);

	return 0;

}
