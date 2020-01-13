#pragma once
# include <vector>
# include "base.h"
#include <set>

Scalar Median(std::vector<Scalar>::iterator begin, std::vector<Scalar>::iterator end);

Scalar compute_points_dist(Point& p1, Point& p2);

void update_maxCorr_and_lag(Scalar max_Dsigma, std::multiset<Scalar>& max_corr, Scalar lag, Point point, std::vector<std::pair<Point, Scalar>>& pairs_point_and_scale);

std::tuple<Point, Point, Point> find_2_farthest_points(Point p0, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& points_profiles);

Point find_point_matching(std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>> & p_profile, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& pts_profiles, VectorType w, Scalar alpha, Scalar nb_samples);

std::vector<std::pair<Point, Point>> compute_matching_pairs_after_rescaling(std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& source_gls_profiles, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& target_gls_profiles, VectorType w, Scalar alpha, Scalar nb_samples, bool cross_check=false);

std::vector<std::tuple<Point, Point, Scalar>> compute_3_closest_pairs_with_scale(std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& source_gls_profiles, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& target_gls_profiles, Scalar ratio, int nb_samples, VectorType w, Scalar alpha = 1.0);

std::vector<std::tuple<Point, Point, Scalar, Scalar>> compute_3_closest_pairs(std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& source_gls_profiles, std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& target_gls_profiles, Scalar ratio, int nb_samples, VectorType w, Scalar alpha = 1.0);

std::pair<Scalar, std::vector<std::tuple<Point, Point, Scalar>>> relative_scale_estimation(std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& source_gls_profiles, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& target_gls_profiles, Scalar ratio, int nb_samples, VectorType w = VectorType::Ones(), Scalar alpha = 1.0);