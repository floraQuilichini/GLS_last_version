#pragma once
#include <vector>
#include <string>
#include <tuple>
#include "base.h"
#include "matching_prioritization.h"
#include "relative_scale.h"


std::string extract_ext(std::string filename);


bool read_ply_file(std::string filename, std::vector<Point>& pointCloud);


bool read_descriptors_text_file(std::string filename, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& gls_profiles, int& nb_points, int& nb_samples, Scalar& min_scale, Scalar& max_scale, Scalar& base, std::map<Point, std::vector<Scalar>>* geom_var_ptr = nullptr);


bool read_descriptors_text_file(std::string filename, std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& gls_profiles_geom_var, int& nb_points, int& nb_samples, Scalar& min_scale, Scalar& max_scale, Scalar& base);


void write_complete_profiles(int nb_points, int nb_samples, Scalar min_scale, Scalar max_scale, Scalar base, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>>>>& points_gls_descriptors_over_scales, std::string output_filename, std::vector<std::pair<Point, std::vector<Scalar>>>* geom_var_ptr = nullptr);


void write_profiles(int nb_points, int nb_samples, Scalar min_scale, Scalar max_scale, Scalar base, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>>>>& points_gls_descriptors_over_scales, std::string output_filename);


void write_points(int nb_points, int nb_samples, Scalar min_scale, Scalar max_scale, Scalar base, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, VectorType, Scalar, Scalar>>>>& points_gls_descriptors_over_scales, std::string output_filename);


void write_relative_scale_and_matching_points(std::string output_filename, std::pair<Scalar, std::vector<std::tuple<Point, Point, Scalar>>>& relative_scale_and_matching_points);

void write_matching_points(std::string output_filename, std::vector<std::pair<Point, Point>>& matching_points);

void write_matrix_transform(Eigen::Matrix4d& transform, std::string& output_filename);

void write_closest_matching_points(Point& target_point, std::vector<std::tuple<Point, Scalar, bool>>& pairs_source_and_scale, std::string& output_filename, bool new_file = false);

void write_closest_matching_points(const Point& target_point, std::vector<std::tuple<Point, Scalar, bool>>& pairs_source_and_scale, std::string& output_filename, bool new_file = false);

void write_closest_matching_points(Point& target_point, std::vector<std::tuple<int, Point, Scalar, Scalar, bool>>& pairs_source_lag_and_corr, std::string& output_filename, bool new_file = false);

void write_closest_matching_points(pair_priority_queue& queue, PointMap* source_pointMap, PointMap* target_pointMap, std::string& output_filename, bool new_file = false);

void write_tuples(std::string& filename, pair_priority_queue queue, PointMap* source_pointMap, PointMap* target_pointMap, Scalar max_err_scale, Scalar bbox_diag, Scalar lambda);

void write_kpairs(std::string& filename, pair_priority_queue queue, PointMap* source_pointMap, PointMap* target_pointMap, Scalar max_err_scale, int nb_pairs, Scalar bbox_size, Scalar lambda);

void write_ply_file_for_debug(std::string& filename, PointMap* target_pointMap, std::priority_queue<std::tuple<int, int, Scalar, Scalar>, std::vector<std::tuple<int, int, Scalar, Scalar>>, mycomparison> k_pair_queue, int index_first_point);