
#include "stdafx.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <string>
#include "Patate/grenaille.h"
#include "Eigen/Eigen"
#include "IO.h"

using namespace std;
using namespace Grenaille;


std::string extract_ext(std::string filename)
{
	size_t i = filename.rfind('.', filename.length());
	if (i != std::string::npos) {
		return(filename.substr(i + 1, filename.length() - i));
	}

	return("");
}


// Define related structure
typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightFunc;
typedef Basket<Point, WeightFunc, UnorientedSphereFit, GLSParam> Fit2;
typedef Basket<Point, WeightFunc, OrientedSphereFit, GLSParam> Fit1;
bool read_ply_file(std::string filename, vector<Point>& pointCloud)
{
	std::string file_ext = extract_ext(filename);
	if (file_ext.compare("ply") == 0)
	{
		// read ply file and go to the section where faces are listed
		std::ifstream file;
		file.open(filename);

		if (!file)
		{
			std::cerr << "File " << filename << " could not be opened" << std::endl;
			return false;
		}

		// find number of vertices
		std::string header_line;
		int nb_vertices = 0;

		do {
			getline(file, header_line);
			//std::cout << header_line << std::endl;
			if (header_line.find("element vertex") != std::string::npos)
				nb_vertices = std::atoi((header_line.erase(0, 15)).c_str());

		} while (header_line.compare("end_header") != 0);

		// get points
		Scalar x, y, z, nx, ny, nz;
		int counter = 0;
		while (counter < nb_vertices)
		{
			file >> x >> y >> z >> nx >> ny >> nz;

			// normalize normals
			Scalar n_norm = sqrt(nx*nx + ny*ny + nz*nz);
			nx = nx / n_norm;
			ny = ny / n_norm;
			nz = nz / n_norm;

			// push point
			pointCloud.push_back(Point({ x, y, z }, { nx, ny, nz }));
			counter = counter + 1;
			//std::cout << "x : " << x << " y : " << y << " z : " << z << std::endl;
		}

		file.close();
		return true;
	}

	else
	{
		// not implemented yet for other types of files
		return false;
	}
}



bool read_descriptors_text_file(std::string filename, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& gls_profiles, int& nb_points, int& nb_samples, Scalar& min_scale, Scalar& max_scale, Scalar& base, std::map<Point, std::vector<Scalar>, PointComp>* geom_var_ptr)
{
	std::string file_ext = extract_ext(filename);
	if (file_ext.compare("txt") == 0)
	{
		// read text file
		std::ifstream file;
		file.open(filename);

		if (!file)
		{
			std::cerr << "File " << filename << " could not be opened" << std::endl;
			return false;
		}

		// get header
		std::string header_line;
		getline(file, header_line);
		min_scale = std::atof(header_line.c_str());
		getline(file, header_line);
		max_scale = std::atof(header_line.c_str());
		getline(file, header_line);
		base = std::atof(header_line.c_str());
		getline(file, header_line);
		nb_samples = std::atoi(header_line.c_str());
		getline(file, header_line);
		nb_points = std::atoi(header_line.c_str());

		// get points and descriptors
		Scalar x, y, z, nx, ny, nz;
		Scalar tau, kappa, phi, geom_var;
		int counter1 = 0;
		
		if (geom_var_ptr)
		{ 
			while (counter1 < nb_points)
			{
				file >> x >> y >> z >> nx >> ny >> nz;
				//std::cout << "point : " << x << " " << y  << " " << z << " " << nx << " " << ny << " " << nz << " " << std::endl;

				int counter2 = 0;
				std::vector<std::tuple<Scalar, Scalar, Scalar>> scale_profiles;
				std::vector<Scalar> geom_variations;
				while (counter2 < nb_samples)
				{
					file >> tau >> kappa >> phi >> geom_var;
					//std::cout << "parameters : " << tau << kappa << phi << std::endl;
					scale_profiles.push_back(std::make_tuple(tau, kappa, phi));
					geom_variations.push_back(geom_var);
					counter2 = counter2 + 1;
				}

				gls_profiles.push_back(std::make_pair(Point({ x, y, z }, { nx, ny, nz }), scale_profiles));
				geom_var_ptr->insert(std::make_pair(Point({ x, y, z }, { nx, ny, nz }), geom_variations));
				counter1 = counter1 + 1;
			}
		}
		else
		{
			while (counter1 < nb_points)
			{
				file >> x >> y >> z >> nx >> ny >> nz;
				//std::cout << "point : " << x << " " << y  << " " << z << " " << nx << " " << ny << " " << nz << " " << std::endl;

				int counter2 = 0;
				std::vector<std::tuple<Scalar, Scalar, Scalar>> scale_profiles;
				while (counter2 < nb_samples)
				{
					file >> tau >> kappa >> phi;
					//std::cout << "parameters : " << tau << kappa << phi << std::endl;
					scale_profiles.push_back(std::make_tuple(tau, kappa, phi));
					counter2 = counter2 + 1;
				}

				gls_profiles.push_back(std::make_pair(Point({ x, y, z }, { nx, ny, nz }), scale_profiles));
				counter1 = counter1 + 1;
			}
		}

		file.close();
		return true;
	}

	else
	{
		// not implemented yet for other types of files
		return false;
	}
}








void write_complete_profiles(int nb_points, int nb_samples, Scalar min_scale, Scalar max_scale, Scalar base, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>>>>& points_gls_descriptors_over_scales, std::string output_filename, std::vector<std::pair<Point, std::vector<Scalar>>>* geom_var_ptr)
{

	std::ofstream output_file;
	output_file.open(output_filename, ios::out);
	// header
	output_file << min_scale << std::endl;
	output_file << max_scale << std::endl;
	output_file << base << std::endl;
	output_file << nb_samples << std::endl;
	output_file << nb_points << std::endl;

	if (geom_var_ptr)
	{
		// create iterator
		std::vector<std::pair<Point, std::vector<Scalar>>>::iterator it = geom_var_ptr->begin();
		// Point-descriptor over scales
		for (int i = 0; i < nb_points; i++)
		{
			Point point = (points_gls_descriptors_over_scales[i]).first;
			output_file << (point.pos()).transpose() << " " << (point.normal()).transpose() << std::endl;
			std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>> descriptors_over_scales = (points_gls_descriptors_over_scales[i]).second;
			for (int j = 0; j < nb_samples; j++)
				output_file << std::get<0>(descriptors_over_scales[j]) << " " << std::get<2>(descriptors_over_scales[j]) << " " << std::get<3>(descriptors_over_scales[j]) <<  " " << it->second[j] <<  std::endl;

			it++;
		}
		
	}
	else
	{
		// Point-descriptor over scales
		for (int i = 0; i < nb_points; i++)
		{
			Point point = (points_gls_descriptors_over_scales[i]).first;
			output_file << (point.pos()).transpose() << " " << (point.normal()).transpose() << std::endl;
			std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>> descriptors_over_scales = (points_gls_descriptors_over_scales[i]).second;
			for (int j = 0; j < nb_samples; j++)
				output_file << std::get<0>(descriptors_over_scales[j]) << " " << std::get<2>(descriptors_over_scales[j]) << " " << std::get<3>(descriptors_over_scales[j]) << std::endl;
		}
	}
	output_file.close();
}



void write_profiles(int nb_points, int nb_samples, Scalar min_scale, Scalar max_scale, Scalar base, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>>>>& points_gls_descriptors_over_scales, std::string output_filename)
{

	std::ofstream output_file;
	output_file.open(output_filename, ios::out);
	// header
	output_file << min_scale << std::endl;
	output_file << max_scale << std::endl;
	output_file << base << std::endl;
	output_file << nb_samples << std::endl;
	output_file << nb_points << std::endl;
	// Point-descriptor over scales
	for (int i = 0; i < nb_points; i++)
	{
		std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>> descriptors_over_scales = (points_gls_descriptors_over_scales[i]).second;
		for (int j = 0; j < nb_samples; j++)
			output_file << std::get<0>(descriptors_over_scales[j]) << " " << std::get<2>(descriptors_over_scales[j]) << " " << std::get<3>(descriptors_over_scales[j]) << std::endl;
	}
	output_file.close();
}



void write_points(int nb_points, int nb_samples, Scalar min_scale, Scalar max_scale, Scalar base, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>>>>& points_gls_descriptors_over_scales, std::string output_filename)
{

	std::ofstream output_file;
	output_file.open(output_filename, ios::out);
	// header
	output_file << min_scale << std::endl;
	output_file << max_scale << std::endl;
	output_file << base << std::endl;
	output_file << nb_samples << std::endl;
	output_file << nb_points << std::endl;
	// Point-descriptor over scales
	for (int i = 0; i < nb_points; i++)
	{
		Point point = (points_gls_descriptors_over_scales[i]).first;
		output_file << (point.pos()).transpose() << " " << (point.normal()).transpose() << std::endl;

	}
	output_file.close();
}



void write_relative_scale_and_matching_points(std::string output_filename, std::pair<Scalar, std::vector<std::tuple<Point, Point, Scalar>>>& relative_scale_and_matching_points)
{
	std::ofstream output_file;
	output_file.open(output_filename, ios::out);
	// relative scale estimation
	output_file << relative_scale_and_matching_points.first << std::endl; 
	// Matching points and related scales
	for (int i = 0; i < relative_scale_and_matching_points.second.size(); i++)
		output_file << std::get<0>(relative_scale_and_matching_points.second[i]).pos().transpose() << " " << std::get<1>(relative_scale_and_matching_points.second[i]).pos().transpose() << " " << std::get<2>(relative_scale_and_matching_points.second[i]) << std::endl;

	output_file.close();
}


void write_matching_points(std::string output_filename, std::vector<std::pair<Point, Point>>& matching_points)
{
	std::ofstream output_file;
	output_file.open(output_filename, ios::out);
	// Matching points
	for (int i = 0; i < matching_points.size(); i++)
		output_file << (matching_points[i].first).pos().transpose() << " " << (matching_points[i].second).pos().transpose() << std::endl;

	output_file.close();
}

