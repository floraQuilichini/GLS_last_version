#pragma once
# include <vector>
# include "base.h"
#include <set>

class PointMap
{
private:
	std::map<int, std::tuple<Point, Eigen::ArrayX3d, Scalar>> map_index_point_profiles_cost_;

public:
	PointMap(){};
	void set_point_profiles_cost(std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& point_gls_profiles, int nb_samples, int alpha = 1.0);
	std::pair<std::map<int, std::tuple<Point, Eigen::ArrayX3d, Scalar>>::iterator, std::map<int, std::tuple<Point, Eigen::ArrayX3d, Scalar>>::iterator> get_iterator_range();
	std::map<int, std::tuple<Point, Eigen::ArrayX3d, Scalar>>::iterator find_index(int index);
};


Scalar Median(std::vector<Scalar>::iterator begin, std::vector<Scalar>::iterator end);

Eigen::Array<Scalar, 1, 3> get_max(const Eigen::Array<Scalar, 1, 3>& v1, const Eigen::Array<Scalar, 1, 3>& v2);

void update_maxCorr_and_lag(Scalar max_Dsigma, std::multiset<Scalar>& max_corr, Scalar lag, Point point, std::vector<std::tuple<Point, Scalar, bool>>& pairs_point_and_scale, Scalar ratio);

void update_maxCorr_and_lag(Scalar max_Dsigma, std::multiset<Scalar>& max_corr, Scalar lag, int point_index, Point& point, std::vector<std::tuple<int, Point, Scalar, Scalar, bool>>& pairs_point_lag_and_corr, Scalar ratio);

Point find_point_matching(std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>> & p_profile, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& pts_profiles, VectorType w, Scalar alpha, Scalar nb_samples);

std::vector<std::pair<Point, Point>> compute_matching_pairs_after_rescaling(std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& source_gls_profiles, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& target_gls_profiles, VectorType w, Scalar alpha, Scalar nb_samples, bool cross_check=false);

std::vector<std::tuple<Point, Point, Scalar>> compute_3_closest_pairs_with_scale(std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& source_gls_profiles, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& target_gls_profiles, Scalar ratio, int nb_samples, VectorType w, Scalar alpha = 1.0);

std::vector<std::tuple<Point, Point, Scalar, Scalar>> compute_3_closest_pairs(std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& source_gls_profiles, std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& target_gls_profiles, Scalar ratio, int nb_samples, VectorType w, Scalar alpha = 1.0);

std::vector<std::tuple<Point, Point, Scalar, Scalar>> compute_3_closest_pairs(std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& source_gls_profiles, std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& target_gls_profiles, Scalar ratio, int nb_source_samples, int nb_target_samples, VectorType w, Scalar alpha = 1.0);

std::vector<std::tuple<std::pair<int, int>, Scalar, Scalar, Scalar>> compute_symmetric_pairs(PointMap& point_map_source, PointMap& point_map_target, Scalar ratio, int k, VectorType w, Scalar alpha = 1.0, bool cross_check = false);

std::pair<Scalar, Scalar> compute_optimal_shift(Eigen::ArrayX3d& source_profiles, Eigen::ArrayX3d& target_profiles, Eigen::Array<Scalar, 1, 3>& W, Scalar alpha, Scalar ratio, bool reverse = false);

std::pair<Scalar, std::vector<std::tuple<Point, Point, Scalar>>> relative_scale_estimation(std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& source_gls_profiles, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& target_gls_profiles, Scalar ratio, int nb_samples, VectorType w = VectorType::Ones(), Scalar alpha = 1.0);