#pragma once
#include <vector>
#include <string>
#include <tuple>
#include "base.h"
#include "Eigen/Eigen"
#include "matching_prioritization.h"

typedef Point::Scalar Scalar;
typedef Point::VectorType VectorType;
class RansacScheme
{
private:
	std::priority_queue<std::tuple<Point, Point, Scalar, Scalar>, std::vector<std::tuple<Point, Point, Scalar, Scalar>>, mycomparison>* queue_ptr_;
	std::priority_queue<std::tuple<Point, Point, Scalar, Scalar>, std::vector<std::tuple<Point, Point, Scalar, Scalar>>, mycomparison> queue_copy_;

public:
	RansacScheme(pair_priority_queue& pair_prioritized);
	struct triplet {
		std::tuple<Point, Point, Scalar, Scalar> pair1;
		std::tuple<Point, Point, Scalar, Scalar> pair2;
		std::tuple<Point, Point, Scalar, Scalar> pair3;
	};
	triplet pop_triplet();
	Scalar scaleDiff(triplet t);
	Eigen::Matrix4d compute_rigid_transform(triplet);
	Scalar registrationErr(Eigen::Matrix4d transform, std::vector<std::pair<Point, Point>>& pairs_source_target);
	Scalar compute_angle(VectorType v1, VectorType v2);
	Scalar normalErr(Eigen::Matrix4d transform, std::vector<std::pair<Point, Point>>& pairs_source_target);
};