#pragma once
#include <vector>
#include <string>
#include <tuple>
#include "base.h"


std::string extract_ext(std::string filename);


bool read_ply_file(std::string filename, std::vector<Point>& pointCloud);


bool read_descriptors_text_file(std::string filename, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& gls_profiles, int& nb_points, int& nb_samples, Scalar& min_scale, Scalar& max_scale, Scalar& base, std::map<Point, std::vector<Scalar>>* geom_var_ptr = nullptr);


bool read_descriptors_text_file(std::string filename, std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& gls_profiles_geom_var, int& nb_points, int& nb_samples, Scalar& min_scale, Scalar& max_scale, Scalar& base);


void write_complete_profiles(int nb_points, int nb_samples, Scalar min_scale, Scalar max_scale, Scalar base, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>>>>& points_gls_descriptors_over_scales, std::string output_filename, std::vector<std::pair<Point, std::vector<Scalar>>>* geom_var_ptr = nullptr);


void write_profiles(int nb_points, int nb_samples, Scalar min_scale, Scalar max_scale, Scalar base, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>>>>& points_gls_descriptors_over_scales, std::string output_filename);


void write_points(int nb_points, int nb_samples, Scalar min_scale, Scalar max_scale, Scalar base, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>>>>& points_gls_descriptors_over_scales, std::string output_filename);


void write_relative_scale_and_matching_points(std::string output_filename, std::pair<Scalar, std::vector<std::tuple<Point, Point, Scalar>>>& relative_scale_and_matching_points);

void write_matching_points(std::string output_filename, std::vector<std::pair<Point, Point>>& matching_points);