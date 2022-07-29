#include "stdafx.h"
#include "math_approx.h"

 tanh_lookup_table::tanh_lookup_table(Scalar min_abscisse, Scalar max_abscisse, int nb_samples)
{
	 min_abscisse_ = min_abscisse;
	 max_abscisse_ = max_abscisse;
	 nb_samples_ = nb_samples;

	step_ = (max_abscisse - min_abscisse) / (nb_samples - 1);

	lookup_table_.reserve(nb_samples);
	for (int i = 0; i < nb_samples; i++)
	{
		Scalar current_abscisse = min_abscisse + i*step_;
		Scalar tanh_value = tanh(current_abscisse);
		lookup_table_[i] = tanh_value;
	}
}


 void tanh_lookup_table::read_lookup_table(const Eigen::ArrayXd& abscisses, Eigen::ArrayXd& tanh_values)
 {
	 size_t size_vec = abscisses.rows();
	 tanh_values = Eigen::ArrayXd(size_vec, 1);

	 // convert abscisse values to int
	 for (int i = 0; i < size_vec; i++)
	 {
		 Scalar value = abscisses(0, i);
		 if (value < min_abscisse_)
			 tanh_values(i) = -1;
		 else
		 {
			 if (value > max_abscisse_)
				 tanh_values(i) = 1;
			 else
			 {
				 int corresponding_int = (int)((value - min_abscisse_) / step_);
				 tanh_values(i) = lookup_table_[corresponding_int];
			 }
		 }
	}
 }

 
 std::vector<VectorType> compute_regular_bbox_dims(Eigen::Matrix3Xd& data)
 {
	 VectorType min_pos, max_pos, dims;
	 min_pos = data.block(0, 0, 3, 1);
	 max_pos = min_pos;

	 for (int i = 0; i < data.cols(); i++)
	 {
		 VectorType point = data.block(0, i, 3, 1);
		 // for x axis
		 if (point(0) < min_pos(0))
			 min_pos(0) = point(0);
		 else
		 {
			 if (point(0) > max_pos(0))
				 max_pos(0) = point(0);
		 }

		 // for y axis
		 if (point(1) < min_pos(1))
			 min_pos(1) = point(1);
		 else
		 {
			 if (point(1) > max_pos(1))
				 max_pos(1) = point(1);
		 }

		 // for z axis
		 if (point(2) < min_pos(2))
			 min_pos(2) = point(2);
		 else
		 {
			 if (point(2) > max_pos(2))
				 max_pos(2) = point(2);
		 }
	 }

	 dims = max_pos - min_pos;
	 std::vector<VectorType> bbox_params;
	 bbox_params.push_back(dims);
	 bbox_params.push_back(min_pos);
	 bbox_params.push_back(max_pos);

	 return bbox_params;
 }



 Eigen::Matrix3d compute_pca_bbox_rescaling(PointMap& pointMap)
 {
	 // load data into eigen matrices
	 Eigen::Matrix3Xd data(3, pointMap.get_size());
	 auto pair_it = pointMap.get_iterator_range();
	 int counter = 0;
	 for (auto it = pair_it.first; it != pair_it.second; it++)
	 {
		 VectorType& pos = std::get<0>(it->second).pos();
		 data.col(counter) = pos;
		 counter +=1;
	 }

	// center data
	 Eigen::Matrix3Xd data_centered = data.colwise() - data.rowwise().mean();

	 // compute covariance
	 Eigen::Matrix3d cov = data_centered * data_centered.adjoint();

	 // get eigen vectors and eigen values
	 Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(cov);
	 Eigen::Matrix<double, 3, 1> eigen_vals = eig.eigenvalues();
	 Eigen::Matrix<double, 3, 3> eigen_vecs = eig.eigenvectors();


	 // the eigen vectors give us a new basis of R3
		 // now we compute the data coordinates in the new base
	Eigen::Matrix3Xd data_centered_in_eigvec_base = eigen_vecs * data_centered;

	// compute "regular" bbox parameters in new frame 
	std::vector<VectorType> bbox_params = compute_regular_bbox_dims(data_centered_in_eigvec_base);

	// compute rescaling matrix
	Eigen::Matrix3d rescaling_mat = Eigen::Matrix3d::Identity();
	rescaling_mat(0, 0) = bbox_params[0](0);
	rescaling_mat(1, 1) = bbox_params[0](1);
	rescaling_mat(2, 2) = bbox_params[0](2);
	
	//rescaling_mat *= eigen_vecs.inverse();

	return rescaling_mat;

 }
 