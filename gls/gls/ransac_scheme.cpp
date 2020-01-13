#include "stdafx.h"

#include <cmath>
#include <algorithm>
#include <vector>
#include <tuple>
#include <string>
#include <iterator>
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

std::pair<Scalar, Scalar> RansacScheme::scaleDiff(RansacScheme::triplet t)
{
	Scalar avgScale = (std::get<2>(t.pair1) + std::get<2>(t.pair2) + std::get<2>(t.pair3)) / 3.0;
	Scalar varScale = (abs(std::get<2>(t.pair1) - avgScale) + abs(std::get<2>(t.pair2) - avgScale) + abs(std::get<2>(t.pair3) - avgScale)) / 3.0;
	return std::make_pair(varScale, avgScale);
}


/*for that function, we follow the tutorial http://nghiaho.com/?page_id=671 */
Eigen::Matrix4d RansacScheme::compute_rigid_transform(RansacScheme::triplet t)
{
	// points 
	VectorType target_pos1 = std::get<0>(t.pair1).pos();
	VectorType source_pos1 = std::get<1>(t.pair1).pos();
	VectorType target_pos2 = std::get<0>(t.pair2).pos();
	VectorType source_pos2 = std::get<1>(t.pair2).pos();
	VectorType target_pos3 = std::get<0>(t.pair3).pos();
	VectorType source_pos3 = std::get<0>(t.pair3).pos();

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

	std::cout << "transform : " << transform << std::endl;

	return transform;
}


/* not sure this is the proper way to solve registration problem for affine transformation with 4 matching points. 
The method used (by solving a linear system) returns a 4-by-4 matrix (with coefficients for rotation, translation and scaling) but solves the 12 unknown as if they were decorrelated 
In fact, it is not quite the case. Some of the terms of the matrix are related to each other and we have only 10 unknowns [1 theta, 3 rot axis, 3 translation, 3 scaling ](and not 12)*/
Eigen::Matrix4d RansacScheme::compute_rigid_transform(RansacScheme::triplet t, std::tuple<Point, Point, Scalar, Scalar> q)
{
	// 4 pairs of points 
	VectorType target_pos1 = std::get<0>(t.pair1).pos();
	VectorType source_pos1 = std::get<1>(t.pair1).pos();
	VectorType target_pos2 = std::get<0>(t.pair2).pos();
	VectorType source_pos2 = std::get<1>(t.pair2).pos();
	VectorType target_pos3 = std::get<0>(t.pair3).pos();
	VectorType source_pos3 = std::get<1>(t.pair3).pos();
	VectorType target_q_pos = std::get<0>(q).pos();
	VectorType source_q_pos = std::get<1>(q).pos();

	// target and source data
	Eigen::Matrix4d source_data = Eigen::Matrix4d::Ones(4, 4);
	Eigen::Matrix4d target_data = Eigen::Matrix4d::Ones(4, 4);
	source_data.block(0, 0, 3, 1) = Eigen::Map<Eigen::ArrayXd>(source_pos1.data(), 3);
	source_data.block(0, 1, 3, 1) = Eigen::Map<Eigen::ArrayXd>(source_pos2.data(), 3);
	source_data.block(0, 2, 3, 1) = Eigen::Map<Eigen::ArrayXd>(source_pos3.data(), 3);
	source_data.block(0, 3, 3, 1) = Eigen::Map<Eigen::ArrayXd>(source_q_pos.data(), 3);

	target_data.block(0, 0, 3, 1) = Eigen::Map<Eigen::ArrayXd>(target_pos1.data(), 3);
	target_data.block(0, 1, 3, 1) = Eigen::Map<Eigen::ArrayXd>(target_pos2.data(), 3);
	target_data.block(0, 2, 3, 1) = Eigen::Map<Eigen::ArrayXd>(target_pos3.data(), 3);
	target_data.block(0, 3, 3, 1) = Eigen::Map<Eigen::ArrayXd>(target_q_pos.data(), 3);

	// solve linear system 
	Eigen::Matrix4d transform = target_data.colPivHouseholderQr().solve(source_data);

	std::cout << "transform : " << transform << std::endl;

	return transform;
}



Scalar RansacScheme::registrationErr(Eigen::Matrix4d transform, std::vector<std::tuple<Point, Point, Scalar, Scalar>>& pairs_source_target)
{
	Scalar err = 0.0;
	int nb_pairs = pairs_source_target.size();
	Eigen::Matrix3d R = transform.block(0, 0, 3, 3);
	Eigen::Vector3d T = transform.block(0, 3, 3, 1);
	for (int k = 0; k < nb_pairs; k++)
		err += ((Eigen::Map<Eigen::MatrixXd>(std::get<0>(pairs_source_target[k]).pos().data(),3, 1) -R*Eigen::Map<Eigen::MatrixXd>(std::get<1>(pairs_source_target[k]).pos().data(), 3, 1) - T).cwiseAbs()).sum();
	
	std::cout << "reg err : " << err / (Scalar)nb_pairs;
	return err / (Scalar)nb_pairs;
}

Scalar RansacScheme::compute_angle(VectorType v1, VectorType v2)
{
	Scalar dot_product = v1.dot(v2);
	Scalar norm_v1 = v1.norm();
	Scalar norm_v2 = v2.norm();
	return acos(dot_product / (norm_v1*norm_v2)); // en radians
}

Scalar RansacScheme::normalErr(Eigen::Matrix4d transform, std::vector<std::tuple<Point, Point, Scalar, Scalar>>& pairs_source_target)
{
	Scalar err = 0.0;
	int nb_pairs = pairs_source_target.size();
	Eigen::Matrix3d R = transform.block(0, 0, 3, 3);
	for (int k = 0; k < nb_pairs; k++)
		err += compute_angle(Eigen::Map<Eigen::MatrixXd>(std::get<0>(pairs_source_target[k]).normal().data(), 3, 1) ,  R*Eigen::Map<Eigen::MatrixXd>(std::get<1>(pairs_source_target[k]).normal().data(), 3, 1));

	std::cout << "  normal err : " << err / (Scalar)nb_pairs << std::endl;
	return err / (Scalar)nb_pairs;
}

bool RansacScheme::is_q_unique(RansacScheme::triplet t, std::tuple<Point, Point, Scalar, Scalar> q)
{
	if (std::get<0>(t.pair1) == std::get<0>(q))
		return false;
	if (std::get<1>(t.pair1) == std::get<1>(q))
		return false;
	return true;
}


bool RansacScheme::is_valid(std::tuple<Point, Point, Scalar, Scalar> q, RansacScheme::triplet t, std::vector<std::tuple<Point, Point, Scalar, Scalar>>& pairs_source_target, Scalar max_err_reg, Scalar max_err_norm)
{
	if (!is_q_unique(t, q))
		return false;
	Eigen::Matrix4d transform = compute_rigid_transform(t, q);
	Scalar reg_err = registrationErr(transform, pairs_source_target);
	Scalar norm_err = normalErr(transform, pairs_source_target);
	if (!(reg_err < max_err_reg && norm_err < max_err_norm))
		return false;
	else
		return true;

}


void RansacScheme::rescale_data(std::vector<std::tuple<Point, Point, Scalar, Scalar>>& pairs_source_target, Scalar avgScale)
{
	for (int k = 0; k < pairs_source_target.size(); k++)
		std::get<1>(pairs_source_target[k]).pos() *= avgScale;
}


Eigen::Matrix4d RansacScheme::ransac_algorithm(int nb_iterations, Scalar max_err_scale, Scalar max_err_reg, Scalar max_err_norm, std::vector<std::tuple<Point, Point, Scalar, Scalar>>& pairs_source_target)
{
	int counter = 0;
	while (counter < nb_iterations)
	{
		triplet t = pop_triplet();
		Scalar err_scale = scaleDiff(t).first;
		std::cout << "error scale " << err_scale <<std::endl;
		if (err_scale < max_err_scale)
		{
			Scalar avgScale = scaleDiff(t).second;
			rescale_data(pairs_source_target, avgScale);
			Eigen::Matrix4d M = compute_rigid_transform(t);
			
			if (registrationErr(M, pairs_source_target) < max_err_reg && normalErr(M, pairs_source_target) < max_err_norm)
			{
				queue_copy_ = *queue_ptr_;
				while (!queue_copy_.empty())
				{
					std::tuple<Point, Point, Scalar, Scalar> q = queue_copy_.top();
					queue_copy_.pop();
					if (is_valid(q, t, pairs_source_target, max_err_reg, max_err_norm))
						return compute_rigid_transform(t, q);
				}
			}

		}
		counter += 1;
	}
	return Eigen::Matrix4d::Identity();
}