#pragma once
#include <vector>
#include <string>
#include <tuple>
#include "base.h"
#include "Eigen/Eigen"
#include "matching_prioritization.h"
#include "relative_scale.h"


class RansacScheme
{
private:
	std::priority_queue<std::tuple<int, int, Scalar, Scalar>, std::vector<std::tuple<int, int, Scalar, Scalar>>, mycomparison>* queue_ptr_;
	std::priority_queue<std::tuple<int, int, Scalar, Scalar>, std::vector<std::tuple<int, int, Scalar, Scalar>>, mycomparison> queue_copy_;
	PointMap* target_pointMap_ = nullptr;
	PointMap* source_pointMap_ = nullptr;

public:
	RansacScheme(pair_priority_queue& pair_prioritized);
	RansacScheme(pair_priority_queue& pair_prioritized, PointMap* source_pointMap, PointMap* target_pointMap);
	struct triplet {
		std::tuple<int, int, Scalar, Scalar> pair1;
		std::tuple<int, int, Scalar, Scalar> pair2;
		std::tuple<int, int, Scalar, Scalar> pair3;
	};
	Scalar compute_points_dist(Point& p1, Point& p2);
	void find_2_farthest_pairs(triplet& triplet);
	void find_2_farthest_pairs_with_respect_to_geomVar(triplet& triplet, Scalar bbox_diag_ratio);
	std::tuple<int, int, Scalar, Scalar> find_pair_to_maximize_triangle_surface(std::tuple<int, int, Scalar, Scalar>& pair1, std::tuple<int, int, Scalar, Scalar>& pair2);
	std::tuple<int, int, Scalar, Scalar> find_farthest_pair(std::tuple<int, int, Scalar, Scalar>& pair);
	triplet pop_triplet();
	triplet pop_3_farthest_pairs();
	triplet pop_3_farthest_pairs(Scalar bbox_diag_ratio);
	triplet pick_triplet(const int ind1, const int ind2, const int ind3);
	std::pair<Scalar, Scalar> scaleDiff(triplet t);
	Scalar compute_scale(triplet& triplet);
	Eigen::Matrix4d compute_rigid_transform(triplet);
	Eigen::Matrix4d compute_rigid_transform(triplet, std::tuple<int, int, Scalar, Scalar> q);
	Scalar registrationErr(Eigen::Matrix4d transform, RansacScheme::triplet& triplet);
	Scalar RansacScheme::registrationErr(Eigen::Matrix4d transform, RansacScheme::triplet& triplet, std::tuple<int, int, Scalar, Scalar>& q);
	Scalar compute_angle(VectorType v1, VectorType v2);
	Scalar RansacScheme::normalErr(Eigen::Matrix4d transform, RansacScheme::triplet& triplet);
	Scalar RansacScheme::normalErr(Eigen::Matrix4d transform, RansacScheme::triplet& triplet, std::tuple<int, int, Scalar, Scalar>& q);
	bool is_q_unique(triplet t, std::tuple<int, int, Scalar, Scalar> q);
	bool is_valid(std::tuple<int, int, Scalar, Scalar> q, triplet t, Scalar max_err_reg, Scalar max_err_norm);
	Eigen::Matrix4d ransac_algorithm(int nb_iterations, Scalar max_err_scale, Scalar max_err_reg, Scalar max_err_norm);
};