#include "stdafx.h"

#include <cmath>
#include <algorithm>
#include <vector>
#include <tuple>
#include <string>
#include "Patate/grenaille.h"
#include "Eigen/Eigen"
#include "ransac_scheme.h"

RansacScheme::RansacScheme(pair_priority_queue& pair_prioritized) : queue_ptr_(pair_prioritized.get_queue_ptr()) {}
RansacScheme::triplet RansacScheme::pop_triplet()
{
	RansacScheme::triplet triplet;

	if (queue_ptr_->size() > 2)
	{
		triplet.pair1 = queue_ptr_->top();
		queue_ptr_->pop();
		triplet.pair2 = queue_ptr_->top();
		queue_ptr_->pop();
		triplet.pair3 = queue_ptr_->top();
		queue_ptr_->pop();
	}
	return triplet;
}

Scalar RansacScheme::scaleDiff(RansacScheme::triplet t)
{
	Scalar avgScale = (std::get<2>(t.pair1) + std::get<2>(t.pair2) + std::get<2>(t.pair3)) / 3.0;
	Scalar varScale = (abs(std::get<2>(t.pair1) - avgScale) + abs(std::get<2>(t.pair2) - avgScale) + abs(std::get<2>(t.pair3) - avgScale)) / 3.0;
	return varScale;
}


/*for that function, we follow the tutorial http://nghiaho.com/?page_id=671 */
Eigen::Matrix4d RansacScheme::compute_rigid_transform(RansacScheme::triplet t)
{
	// points 
	VectorType source_pos1 = std::get<0>(t.pair1).pos();
	VectorType target_pos1 = std::get<1>(t.pair1).pos();
	VectorType source_pos2 = std::get<0>(t.pair2).pos();
	VectorType target_pos2 = std::get<1>(t.pair2).pos();
	VectorType source_pos3 = std::get<0>(t.pair3).pos();
	VectorType target_pos3 = std::get<0>(t.pair3).pos();

	// target and source data
	Eigen::Matrix3d source_data, target_data;
	source_data.col(0) = Eigen::Map<Eigen::ArrayXd>(source_pos1.data(), 3);
	source_data.col(1) = Eigen::Map<Eigen::ArrayXd>(source_pos2.data(), 3);
	source_data.col(2) = Eigen::Map<Eigen::ArrayXd>(source_pos3.data(), 3);
	target_data.col(0) = Eigen::Map<Eigen::ArrayXd>(target_pos1.data(), 3);
	target_data.col(1) = Eigen::Map<Eigen::ArrayXd>(target_pos2.data(), 3);
	target_data.col(2) = Eigen::Map<Eigen::ArrayXd>(target_pos3.data(), 3);

	//computing centroids of both datasets
	Eigen::Vector3d centroid_source, centroid_target;
	centroid_source << (source_pos1[0] + source_pos2[0] + source_pos3[0]) / 3.0, (source_pos1[1] + source_pos2[1] + source_pos3[1]) / 3.0, (source_pos1[2] + source_pos2[2] + source_pos3[2]) / 3.0;
	centroid_target << (target_pos1[0] + target_pos2[0] + target_pos3[0]) / 3.0, (target_pos1[1] + target_pos2[1] + target_pos3[1]) / 3.0, (target_pos1[2] + target_pos2[2] + target_pos3[2]) / 3.0;

	// shift datasets to their respective centers and find rotation 
	Eigen::Matrix3d H, U, V, R;
	H = (target_data.colwise() - centroid_target)*(source_data.colwise() - centroid_source).transpose();
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
	U = svd.matrixU();
	V = svd.matrixV();
	R = V*U.transpose();
		//handling special reflection case
	if (R.determinant() < 0.0)
	{
		V.col(2) = -1.0*V.col(2);
		R = V*U.transpose();
	}

	// find translation
	Eigen::Vector3d translation;
	translation = source_data.col(0) - R*target_data.col(0); 

	// global transform 
	Eigen::Matrix4d transform = Eigen::Matrix4d::Zero(4, 4);
	transform.block(0, 0, 3, 3) = R;
	transform.block(0, 3, 3, 1) = translation;
	transform(3, 3) = 1.0;

	return transform;
}