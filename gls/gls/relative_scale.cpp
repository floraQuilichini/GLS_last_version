#include "stdafx.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <iterator>
#include <numeric>
#include <set>
#include <ctime>
#include "Eigen/Eigen"
#include "relative_scale.h"
#include "IO.h"
#include "matching_prioritization.h"
#include "math_approx.h"

using namespace std;


Scalar Median(std::vector<Scalar>::iterator begin, std::vector<Scalar>::iterator end)
{
	const auto size = std::distance(begin, end);
	std::nth_element(begin, begin + size / 2, end);
	return *std::next(begin, size / 2);
}


Eigen::Array<Scalar, 1, 3> get_max(const Eigen::Array<Scalar, 1, 3>& v1, const Eigen::Array<Scalar, 1, 3>& v2)
{
	Eigen::Array<Scalar, 1, 3> maxArray;
	for (int i = 0; i < 3; i++)
		maxArray(i) = std::max(v1(i), v2(i));

	return maxArray;
}

void update_maxCorr_and_lag(Scalar max_Dsigma, std::multiset<Scalar>& max_corr, Scalar lag, Point point, std::vector<std::tuple<Point, Scalar, bool>>& pairs_point_and_scale, Scalar ratio)
{

	/*std::cout << "intitial max_corr : " << std::endl;
	for (auto itt = max_corr.begin(); itt != max_corr.end(); itt++)
	{
		std::cout << *itt << std::endl;
	}*/

	std::multiset<Scalar>::iterator it = max_corr.insert(max_Dsigma);
	int dist = std::distance(max_corr.begin(), it);
	max_corr.erase(max_corr.begin());

	if (dist > 0)
	{
		int count = 0;
		while (count < dist -1)
		{
			pairs_point_and_scale[count] = pairs_point_and_scale[count+1];
			count = count + 1;
		}
		if (*std::next(max_corr.begin(), dist - 1) >= ratio)
			pairs_point_and_scale[dist -1] = std::make_tuple(point, lag, true);
		else
			pairs_point_and_scale[dist - 1] = std::make_tuple(point, lag, false);

		/*std::cout << "max dsigma : " << max_Dsigma << std::endl;
		std::cout << "dist : " << dist << std::endl;
		std::cout << "lag : " << lag << std::endl;*/
	}

	/*std::cout << "final max_corr : " << std::endl;
	for (auto itt = max_corr.begin(); itt != max_corr.end(); itt++)
	{
		std::cout << *itt << std::endl;
	}*/
}

void update_maxCorr_and_lag(Scalar max_Dsigma, std::multiset<Scalar>& max_corr, Scalar lag, int point_index, Point& point, std::vector<std::tuple<int, Point, Scalar, Scalar, bool>>& pairs_point_lag_and_corr, Scalar ratio)
{
	std::multiset<Scalar>::iterator it = max_corr.insert(max_Dsigma);
	int dist = std::distance(max_corr.begin(), it);
	max_corr.erase(max_corr.begin());

	if (dist > 0)
	{
		int count = 0;
		while (count < dist - 1)
		{
			pairs_point_lag_and_corr[count] = pairs_point_lag_and_corr[count + 1];
			count = count + 1;
		}
		if (*std::next(max_corr.begin(), dist - 1) >= ratio)
			pairs_point_lag_and_corr[dist - 1] = std::make_tuple(point_index, point, lag, max_Dsigma, true);
		else
			pairs_point_lag_and_corr[dist - 1] = std::make_tuple(point_index, point, lag, max_Dsigma, false);

	}


}





Point find_point_matching(std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>> & p_profile, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& pts_profiles, VectorType w, Scalar alpha, Scalar nb_samples)
{
	std::vector<Scalar> p_tau_profile, p_kappa_profile, p_phi_profile;

	// get  p profile in separate vectors
	std::transform(std::begin(p_profile.second), std::end(p_profile.second), std::back_inserter(p_tau_profile), [](auto const& tuple) { return std::get<0>(tuple); });
	std::transform(std::begin(p_profile.second), std::end(p_profile.second), std::back_inserter(p_kappa_profile), [](auto const& tuple) { return std::get<1>(tuple); });
	std::transform(std::begin(p_profile.second), std::end(p_profile.second), std::back_inserter(p_phi_profile), [](auto const& tuple) { return std::get<2>(tuple); });

	// convert vector values to matrices
	Eigen::ArrayX3d profile = Eigen::ArrayX3d::Zero(nb_samples, 3);
	profile.col(0) = Eigen::Map<Eigen::ArrayXd>(p_tau_profile.data(), nb_samples);
	profile.col(1) = Eigen::Map<Eigen::ArrayXd>(p_kappa_profile.data(), nb_samples);
	profile.col(2) = Eigen::Map<Eigen::ArrayXd>(p_phi_profile.data(), nb_samples);



	Scalar max_corr = 0.0;
	Point matching_pair;
	for (std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>::iterator it = pts_profiles.begin(); it != pts_profiles.end(); ++it)
	{
		std::vector<std::tuple<Scalar, Scalar, Scalar>> pt_profile = it->second;
		std::vector<Scalar> pt_tau_profile, pt_kappa_profile, pt_phi_profile;

		// get all gls parameter profiles in separate vectors
		std::transform(std::begin(pt_profile), std::end(pt_profile), std::back_inserter(pt_tau_profile), [](auto const& tuple) { return std::get<0>(tuple); });
		std::transform(std::begin(pt_profile), std::end(pt_profile), std::back_inserter(pt_kappa_profile), [](auto const& tuple) { return std::get<1>(tuple); });
		std::transform(std::begin(pt_profile), std::end(pt_profile), std::back_inserter(pt_phi_profile), [](auto const& tuple) { return std::get<2>(tuple); });

		// convert vector values to matrices
		Eigen::ArrayX3d pt_profiles = Eigen::ArrayX3d::Zero(nb_samples, 3);
		pt_profiles.col(0) = Eigen::Map<Eigen::ArrayXd>(pt_tau_profile.data(), nb_samples);
		pt_profiles.col(1) = Eigen::Map<Eigen::ArrayXd>(pt_kappa_profile.data(), nb_samples);
		pt_profiles.col(2) = Eigen::Map<Eigen::ArrayXd>(pt_phi_profile.data(), nb_samples);

		// compute Dsigma
		Eigen::ArrayX3d diff = Eigen::ArrayX3d::Zero(nb_samples, 3);

		diff.block(0, 0, nb_samples, 3) = profile.block(0, 0, nb_samples, 3) - pt_profiles.block(0, 0, nb_samples, 3);
		Eigen::ArrayX3d diff2 = (diff*diff).rowwise()*Eigen::Array3d(w[0], w[1], w[2]).transpose();
		Eigen::ArrayXd diss_over_scales = (diff2).rowwise().sum();
		Eigen::ArrayXd sigma = 1.0 - (alpha*diss_over_scales).tanh();
		Scalar Dsigma = 1.0 / (Scalar)(nb_samples - 1) * sigma.sum();

		// compare max_Dsigma to max_corr to know if ps and pt are matching pair
		if (Dsigma > max_corr)
		{
			max_corr = Dsigma;
			matching_pair = it->first;
		}
	}

	return matching_pair;
}



std::vector<std::pair<Point, Point>> compute_matching_pairs_after_rescaling(std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& source_gls_profiles, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& target_gls_profiles, VectorType w, Scalar alpha, Scalar nb_samples, bool crosscheck)
{
	size_t nb_source_points = source_gls_profiles.size();
	size_t nb_target_points = target_gls_profiles.size();

	if (nb_source_points >= nb_target_points)   // if nb_source_points >= nb_target_points, then match target points on source points
	{
		std::vector<std::pair<Point, Point>> target_source_matchings;

		// match target points on source points
		for (std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>::iterator it_t = target_gls_profiles.begin(); it_t != target_gls_profiles.end(); ++it_t)
		{
			std::vector<std::tuple<Scalar, Scalar, Scalar>> pt_target_profile = it_t->second;
			std::vector<Scalar> pt_tau_profile, pt_kappa_profile, pt_phi_profile;

			// get all gls parameter profiles in separate vectors
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_tau_profile), [](auto const& tuple) { return std::get<0>(tuple); });
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_kappa_profile), [](auto const& tuple) { return std::get<1>(tuple); });
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_phi_profile), [](auto const& tuple) { return std::get<2>(tuple); });

			// convert vector values to matrices
			Eigen::ArrayX3d pt_profiles = Eigen::ArrayX3d::Zero(nb_samples, 3);
			pt_profiles.col(0) = Eigen::Map<Eigen::ArrayXd>(pt_tau_profile.data(), nb_samples);
			pt_profiles.col(1) = Eigen::Map<Eigen::ArrayXd>(pt_kappa_profile.data(), nb_samples);
			pt_profiles.col(2) = Eigen::Map<Eigen::ArrayXd>(pt_phi_profile.data(), nb_samples);

			// find, for each target point its closest point in source points
			Point source_pair;
			std::vector<std::tuple<Scalar, Scalar, Scalar>> ps_matching_profiles;
			Scalar max_corr = 0.0;

			for (std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>::iterator it_s = source_gls_profiles.begin(); it_s != source_gls_profiles.end(); ++it_s)
			{
				std::vector<std::tuple<Scalar, Scalar, Scalar>> ps_source_profile = it_s->second;
				std::vector<Scalar> ps_tau_profile, ps_kappa_profile, ps_phi_profile;

				// get all gls parameter profiles in separate vectors
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_tau_profile), [](auto const& tuple) { return std::get<0>(tuple); });
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_kappa_profile), [](auto const& tuple) { return std::get<1>(tuple); });
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_phi_profile), [](auto const& tuple) { return std::get<2>(tuple); });

				// convert vector values to matrices
				Eigen::ArrayX3d ps_profiles = Eigen::ArrayX3d::Zero(nb_samples, 3);
				ps_profiles.col(0) = Eigen::Map<Eigen::ArrayXd>(ps_tau_profile.data(), nb_samples);
				ps_profiles.col(1) = Eigen::Map<Eigen::ArrayXd>(ps_kappa_profile.data(), nb_samples);
				ps_profiles.col(2) = Eigen::Map<Eigen::ArrayXd>(ps_phi_profile.data(), nb_samples);

				// compute Dsigma
				Eigen::ArrayX3d diff = Eigen::ArrayX3d::Zero(nb_samples, 3);

				diff.block(0, 0, nb_samples, 3) = ps_profiles.block(0, 0, nb_samples, 3) - pt_profiles.block(0, 0, nb_samples, 3);
				Eigen::ArrayX3d diff2 = (diff*diff).rowwise()*Eigen::Array3d(w[0], w[1], w[2]).transpose();
				Eigen::ArrayXd diss_over_scales = (diff2).rowwise().sum();
				Eigen::ArrayXd sigma = 1.0 - (alpha*diss_over_scales).tanh();
				Scalar Dsigma = 1.0 / (Scalar)(nb_samples - 1) * sigma.sum();

				// compare max_Dsigma to max_corr to know if ps and pt are matching pair
				if (Dsigma > max_corr)
				{
					max_corr = Dsigma;
					source_pair = it_s->first;
					ps_matching_profiles = ps_source_profile;
				}
			}

			if (!crosscheck)
				// store matching target-source matching pairs
				target_source_matchings.push_back(std::make_pair(it_t->first, source_pair));
			else
			{
				// perform cross-check (check source->target matching)
				Point target_matching = find_point_matching(std::make_pair(source_pair, ps_matching_profiles), target_gls_profiles, w, alpha, nb_samples);

				if (target_matching.pos() == (it_t->first).pos())
				{	// store matching target-source matching pairs
					target_source_matchings.push_back(std::make_pair(it_t->first, source_pair));
				}
			}
		}
		return target_source_matchings;
	}
	else
		return compute_matching_pairs_after_rescaling(target_gls_profiles, source_gls_profiles, w, alpha, nb_samples);
}


std::vector<std::tuple<Point, Point, Scalar>> compute_3_closest_pairs_with_scale(std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& source_gls_profiles, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& target_gls_profiles, Scalar ratio, int nb_samples, VectorType w, Scalar alpha)
{

	size_t nb_source_points = source_gls_profiles.size();
	size_t nb_target_points = target_gls_profiles.size();
	Scalar ratio_corr = 0.7;

	if (nb_source_points >= nb_target_points)   // if nb_source_points >= nb_target_points, then match target points on source points
	{
		std::vector<std::tuple<Point, Point, Scalar>> target_source_scale_matchings;

		for (std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>::iterator it_t = target_gls_profiles.begin(); it_t != target_gls_profiles.end(); ++it_t)
		{
			std::vector<std::tuple<Scalar, Scalar, Scalar>> pt_target_profile = it_t->second;
			std::vector<Scalar> pt_tau_profile, pt_kappa_profile, pt_phi_profile;

			// get all gls parameter profiles in separate vectors
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_tau_profile), [](auto const& tuple) { return std::get<0>(tuple); });
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_kappa_profile), [](auto const& tuple) { return std::get<1>(tuple); });
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_phi_profile), [](auto const& tuple) { return std::get<2>(tuple); });

			// convert vector values to matrices
			Eigen::ArrayX3d pt_profiles = Eigen::ArrayX3d::Zero(nb_samples, 3);
			pt_profiles.col(0) = Eigen::Map<Eigen::ArrayXd>(pt_tau_profile.data(), nb_samples);
			pt_profiles.col(1) = Eigen::Map<Eigen::ArrayXd>(pt_kappa_profile.data(), nb_samples);
			pt_profiles.col(2) = Eigen::Map<Eigen::ArrayXd>(pt_phi_profile.data(), nb_samples);

			// find, for each target point its closest point in source points and its related estimated scale
			std::vector<std::tuple<Point, Scalar, bool>> pairs_source_and_scale{ std::make_tuple(Point(), 0.0, false),  std::make_tuple(Point(), 0.0, false) , std::make_tuple(Point(), 0.0, false) };
			std::multiset<Scalar> max_corr{ 0.0, 0.0, 0.0 };

			for (std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>::iterator it_s = source_gls_profiles.begin(); it_s != source_gls_profiles.end(); ++it_s)
			{
				std::vector<std::tuple<Scalar, Scalar, Scalar>> ps_source_profile = it_s->second;
				std::vector<Scalar> ps_tau_profile, ps_kappa_profile, ps_phi_profile;

				// get all gls parameter profiles in separate vectors
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_tau_profile), [](auto const& tuple) { return std::get<0>(tuple); });
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_kappa_profile), [](auto const& tuple) { return std::get<1>(tuple); });
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_phi_profile), [](auto const& tuple) { return std::get<2>(tuple); });

				// convert vector values to matrices
				Eigen::ArrayX3d ps_profiles = Eigen::ArrayX3d::Zero(nb_samples, 3);
				ps_profiles.col(0) = Eigen::Map<Eigen::ArrayXd>(ps_tau_profile.data(), nb_samples);
				ps_profiles.col(1) = Eigen::Map<Eigen::ArrayXd>(ps_kappa_profile.data(), nb_samples);
				ps_profiles.col(2) = Eigen::Map<Eigen::ArrayXd>(ps_phi_profile.data(), nb_samples);


				// find relative scale 
				Scalar max_Dsigma =  0.0;
				Scalar lag = 0.0;
				// compute shifts
					Scalar max_shift = nb_samples - floor(nb_samples*ratio);
					// right shifts
					for (int shift = 0; shift < max_shift - 2; shift++)
					{
						Eigen::ArrayX3d diff = Eigen::ArrayX3d::Zero(nb_samples - shift, 3);
						diff.block(0, 0, nb_samples - shift, 3) = ps_profiles.block(0, 0, nb_samples - shift, 3) - pt_profiles.block(shift, 0, nb_samples, 3);
						Eigen::ArrayX3d diff2 = (diff*diff).rowwise()*Eigen::Array3d(w[0], w[1], w[2]).transpose();
						Eigen::ArrayXd diss_over_scales = (diff2).rowwise().sum();
						Eigen::ArrayXd sigma = 1.0 - (alpha*diss_over_scales).tanh();
						Scalar Dsigma = 1.0 / (Scalar)(nb_samples - 1 - shift) * sigma.sum();

						if (Dsigma > max_Dsigma)
						{
							max_Dsigma = Dsigma;
							lag = -shift;
						}

					}
					// left shifts
					for (int shift = 0; shift < max_shift - 2; shift++)
					{
						Eigen::ArrayX3d diff = Eigen::ArrayX3d::Zero(nb_samples - shift, 3);

						diff.block(0, 0, nb_samples - shift, 3) = ps_profiles.block(shift, 0, nb_samples, 3) - pt_profiles.block(0, 0, nb_samples - shift, 3);
						Eigen::ArrayX3d diff2 = (diff*diff).rowwise()*Eigen::Array3d(w[0], w[1], w[2]).transpose();
						Eigen::ArrayXd diss_over_scales = (diff2).rowwise().sum();
						Eigen::ArrayXd sigma = 1.0 - (alpha*diss_over_scales).tanh();
						Scalar Dsigma = 1.0 / (Scalar)(nb_samples - 1 - shift) * sigma.sum();

						if (Dsigma > max_Dsigma)
						{
							max_Dsigma = Dsigma;
							lag = shift;
						}
				
					}

				// compare max_Dsigma to max_corr to know if ps and pt belong to one of the 3 matching pairs
				update_maxCorr_and_lag(max_Dsigma, max_corr, lag, it_s->first, pairs_source_and_scale, ratio_corr);

			}

			// store matching target-source matching pairs
			for (int k =0; k<3; k++ )
				target_source_scale_matchings.push_back(std::make_tuple(it_t->first, std::get<0>(pairs_source_and_scale[k]), std::get<1>(pairs_source_and_scale[k])));
		}

		return target_source_scale_matchings;
	}

	else // otherwise match source points on target points
	{
		return compute_3_closest_pairs_with_scale(target_gls_profiles, source_gls_profiles, ratio, nb_samples, w, alpha);
	}

}


std::vector<std::tuple<Point, Point, Scalar, Scalar>> compute_3_closest_pairs(std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& source_gls_profiles, std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& target_gls_profiles, Scalar ratio, int nb_samples, VectorType w, Scalar alpha)
{
	size_t nb_source_points = source_gls_profiles.size();
	size_t nb_target_points = target_gls_profiles.size();
	Scalar ratio_corr = 0.7;

	if (nb_source_points >= nb_target_points)   // if nb_source_points >= nb_target_points, then match target points on source points
	{
		std::vector<std::tuple<Point, Point, Scalar, Scalar>> target_source_matchings;

		//compute point priority
		std::map<Point, Scalar> target_priority, source_priority;
//#pragma omp parallel for
		for (int k=0; k<nb_target_points; k++)
		{
			Scalar priority = compute_point_priority(std::get<2>(target_gls_profiles[k]), nb_samples, alpha);
			target_priority.insert(std::make_pair(std::get<0>(target_gls_profiles[k]), priority));
		}

		//std::cout << "source cost : " << std::endl;
//#pragma omp parallel for
		for (int k = 0; k<nb_source_points; k++)
		{
			Scalar priority = compute_point_priority(std::get<2>(source_gls_profiles[k]), nb_samples, alpha);
			source_priority.insert(std::make_pair(std::get<0>(source_gls_profiles[k]), priority));
			//std::cout << priority << std::endl;
		}

		// build lookup table
		tanh_lookup_table lookup_table(-2.0, 2.0, 400);

		//max shift value
		int max_shift = nb_samples - floor(nb_samples*ratio);

		// get W vector   -- slower than [w(0), w(1), w(2)].transpose() ?
		Eigen::Array<Scalar,1,3> W;
		W  << w[0], w[1], w[2];
		/*std::cout << "size W : " << W.rows() << " rows and : " << W.cols() << " cols" << std::endl;*/

		// search for target->source matching
//#pragma omp parallel for
		for (int u = 0 ; u < target_gls_profiles.size() ; ++u)
		// for (std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>::iterator it_t = target_gls_profiles.begin(); it_t != target_gls_profiles.end(); ++it_t)
		{
			//std::cout << " p_target = " << std::get<0>(target_gls_profiles[u]).pos().transpose() << std::endl;
			const std::vector<std::tuple<Scalar, Scalar, Scalar>> & pt_target_profile = std::get<1>(target_gls_profiles[u]);
			// std::vector<std::tuple<Scalar, Scalar, Scalar>> pt_target_profile = std::get<1>(*it_t);
			std::vector<Scalar> pt_tau_profile, pt_kappa_profile, pt_phi_profile;

			// get all gls parameter profiles in separate vectors
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_tau_profile), [](auto const& tuple) { return std::get<0>(tuple); });
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_kappa_profile), [](auto const& tuple) { return std::get<1>(tuple); });
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_phi_profile), [](auto const& tuple) { return std::get<2>(tuple); });

			// convert vector values to matrices
			Eigen::ArrayX3d pt_profiles(nb_samples, 3);
			pt_profiles.col(0) = Eigen::Map<Eigen::ArrayXd>(pt_tau_profile.data(), nb_samples);
			pt_profiles.col(1) = Eigen::Map<Eigen::ArrayXd>(pt_kappa_profile.data(), nb_samples);
			pt_profiles.col(2) = Eigen::Map<Eigen::ArrayXd>(pt_phi_profile.data(), nb_samples);

			// find, for each target point its closest point in source points and its related estimated scale
			std::vector<std::tuple<Point, Scalar, bool>> pairs_source_and_scale{ std::make_tuple(Point(), 0.0, false),  std::make_tuple(Point(), 0.0, false) , std::make_tuple(Point(), 0.0, false) };
			std::multiset<Scalar> max_corr{ 0.0, 0.0, 0.0 };
			//clock_t tstart, tend;
			//tstart = clock();
//#pragma omp parallel for
			for (int v = 0; v < source_gls_profiles.size(); ++v)
			{
				//std::cout << " p_source = " << std::get<0>(source_gls_profiles[v]).pos().transpose() << std::endl;
				const std::vector<std::tuple<Scalar, Scalar, Scalar>> & ps_source_profile = std::get<1>(source_gls_profiles[v]);
				std::vector<Scalar> ps_tau_profile, ps_kappa_profile, ps_phi_profile;

				// get all gls parameter profiles in separate vectors
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_tau_profile), [](auto const& tuple) { return std::get<0>(tuple); });
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_kappa_profile), [](auto const& tuple) { return std::get<1>(tuple); });
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_phi_profile), [](auto const& tuple) { return std::get<2>(tuple); });

				// convert vector values to matrices
				Eigen::ArrayX3d ps_profiles(nb_samples, 3);
				ps_profiles.col(0) = Eigen::Map<Eigen::ArrayXd>(ps_tau_profile.data(), nb_samples);
				ps_profiles.col(1) = Eigen::Map<Eigen::ArrayXd>(ps_kappa_profile.data(), nb_samples);
				ps_profiles.col(2) = Eigen::Map<Eigen::ArrayXd>(ps_phi_profile.data(), nb_samples);


				// find relative scale 
				Scalar max_Dsigma = 0.0;
				Scalar lag = 0.0;
				Eigen::ArrayXd tanh_vec;

				// compute shifts
					// left shifts (positive ones)
//#pragma omp parallel for
					for (int shift = 0; shift < max_shift - 2; shift++)
					{
						Eigen::ArrayX3d diff(nb_samples - shift, 3);

						diff.block(0, 0, nb_samples - shift, 3) = ps_profiles.block(shift, 0, nb_samples - shift, 3) - pt_profiles.block(0, 0, nb_samples - shift, 3);
						Eigen::Array<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> diff2 = (diff*diff).rowwise()*W;
						Eigen::Array<Scalar, Eigen::Dynamic, 1, 1, Eigen::RowMajor> diss_over_scales = (diff2).rowwise().sum();
						Eigen::ArrayXd sigma = 1.0 - (alpha*diss_over_scales).tanh();
						/*lookup_table.read_lookup_table(alpha*diss_over_scales, tanh_vec);
						Eigen::ArrayXd sigma = 1.0 - tanh_vec;*/
						Scalar Dsigma = 1.0 / (Scalar)(nb_samples - 1 - shift) * sigma.sum();

						if (Dsigma > max_Dsigma)
						{
							max_Dsigma = Dsigma;
							lag = shift;
							//std::cout << " max Dsigma : " << max_Dsigma << std::endl;
						}

						/*if (std::get<0>(source_gls_profiles[v]).pos() == std::get<0>(target_gls_profiles[u]).pos())
						{
						std::cout << " Dsigma = " << Dsigma << std::endl;
						std::cout << v << std::endl;
						}*/

						/*if (v == 283 && shift == 25)
						{
							std::cout << "Diff : " << diff << std::endl;
							std::cout << "Diff2 : " << diff2 << std::endl;
							std::cout << " Dsigma : " << Dsigma << std::endl;
						}*/

					}


					// right shifts (negative ones)
//#pragma omp parallel for
				for (int shift = 0; shift < max_shift - 2; shift++)
				{
					Eigen::ArrayX3d diff(nb_samples - shift, 3);
					//std::cout << "diff nb rows : " << diff.rows() << " diff nb cols : " << diff.cols() << std::endl;
					diff.block(0, 0, nb_samples - shift, 3) = ps_profiles.block(0, 0, nb_samples - shift, 3) - pt_profiles.block(shift, 0, nb_samples - shift, 3);
					//std::cout << " diff size : " << diff.size() << std::endl;
					Eigen::Array<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> diff2 = (diff*diff).rowwise()*W;
					//std::cout << " diff2 : " << diff2.size() << std::endl;
					Eigen::Array<Scalar, Eigen::Dynamic, 1, 1, Eigen::RowMajor> diss_over_scales = (diff2).rowwise().sum();
					//std::cout << " diss_over_scales : " << diss_over_scales.size() << std::endl;
					Eigen::ArrayXd sigma = 1.0 - (alpha*diss_over_scales).tanh();
					/*lookup_table.read_lookup_table(alpha*diss_over_scales, tanh_vec);  // look-up table is faster than calling tanh()
					Eigen::ArrayXd sigma = 1.0 - tanh_vec;*/
					//std::cout << " sigma : " << sigma.size() << std::endl;
					Scalar Dsigma = 1.0 / (Scalar)(nb_samples - 1 - shift) * sigma.sum();

					if (Dsigma > max_Dsigma)
					{
						max_Dsigma = Dsigma;
						lag = -shift;
						//std::cout << " max Dsigma : " << max_Dsigma << std::endl;
					}

					/*if (std::get<0>(source_gls_profiles[v]) == std::get<0>(target_gls_profiles[u]))
					{
						std::cout << " Dsigma = " << Dsigma << std::endl;
						std::cout << v << std::endl;
					}*/
					/*if ( v == 283 && shift == 25)
					{
						std::cout << "Diff : " << diff << std::endl;
						std::cout << "Diff2 : " << diff2 << std::endl;
						std::cout << " Dsigma : " << Dsigma << std::endl;
					}*/

				}

				// compare max_Dsigma to max_corr to know if ps and pt belong to one of the 3 matching pairs
				update_maxCorr_and_lag(max_Dsigma, max_corr, lag, std::get<0>(source_gls_profiles[v]), pairs_source_and_scale, ratio_corr);

			}

		/*	// display pair_source_and_scale
			for (int k = 0; k < pairs_source_and_scale.size(); k++)
				std::cout << "source matching point : " << pairs_source_and_scale[k].first.pos() << ", lag : " << pairs_source_and_scale[k].second << std::endl;
		*/	
			std::string output_filename = "C:\\Registration\\test_gls_algo\\matching_pairs\\source_matching_pairs_unsorted.txt";
			write_closest_matching_points(std::get<0>(target_gls_profiles[u]), pairs_source_and_scale, output_filename, false);
		

			// store matching target-source matching pairs with their cost
			for (int k = 0; k < 3; k++)
			{
				const Point & pt = std::get<0>(target_gls_profiles[u]);
				Scalar cost_pair = target_priority.find(pt)->second * source_priority.find(std::get<0>(pairs_source_and_scale[k]))->second;
				//std::cout << "target_cost : " << target_priority.find(pt)->second << ", source_cost : "<< source_priority.find(pairs_source_and_scale[k].first)->second << ", pair cost : "<< cost_pair << std::endl;
				target_source_matchings.push_back(std::make_tuple(pt, std::get<0>(pairs_source_and_scale[k]), std::get<1>(pairs_source_and_scale[k]), cost_pair));
			}

			//tend = clock()-tstart;
			//cout << "time for computing 3 matching pair : " << (float)tend/CLOCKS_PER_SEC << " second(s)." << endl;
		}

		return target_source_matchings;
	}

	else // otherwise match source points on target points
	{
		return compute_3_closest_pairs(target_gls_profiles, source_gls_profiles, ratio, nb_samples, w, alpha);
	}
}



std::vector<std::tuple<Point, Point, Scalar, Scalar>> compute_3_closest_pairs(std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& source_gls_profiles, std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& target_gls_profiles, Scalar ratio, int nb_source_samples, int nb_target_samples, VectorType w, Scalar alpha)
{
	size_t nb_source_points = source_gls_profiles.size();
	size_t nb_target_points = target_gls_profiles.size();
	Scalar ratio_corr = 0.9;

	if (nb_source_points >= nb_target_points)   // if nb_source_points >= nb_target_points, then match target points on source points
	{
		std::vector<std::tuple<Point, Point, Scalar, Scalar>> target_source_matchings;

		//compute point priority
		std::map<Point, Scalar> target_priority, source_priority;
//#pragma omp parallel for
		for (int k = 0; k<nb_target_points; k++)
		{
			Scalar priority = compute_point_priority(std::get<2>(target_gls_profiles[k]), nb_target_samples, alpha);
			target_priority.insert(std::make_pair(std::get<0>(target_gls_profiles[k]), priority));
		}

		//std::cout << "source cost : " << std::endl;
//#pragma omp parallel for
		for (int k = 0; k<nb_source_points; k++)
		{
			Scalar priority = compute_point_priority(std::get<2>(source_gls_profiles[k]), nb_source_samples, alpha);
			source_priority.insert(std::make_pair(std::get<0>(source_gls_profiles[k]), priority));
		}

		// get W vector   -- slower than [w(0), w(1), w(2)].transpose() ?
		Eigen::Array<Scalar, 1, 3> W;
		W << w[0], w[1], w[2];

		// search for target->source matching
//#pragma omp parallel for
		for (int u = 0; u < target_gls_profiles.size(); ++u)
		{
			const std::vector<std::tuple<Scalar, Scalar, Scalar>> & pt_target_profile = std::get<1>(target_gls_profiles[u]);
			std::vector<Scalar> pt_tau_profile, pt_kappa_profile, pt_phi_profile;

			// get all gls parameter profiles in separate vectors
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_tau_profile), [](auto const& tuple) { return std::get<0>(tuple); });
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_kappa_profile), [](auto const& tuple) { return std::get<1>(tuple); });
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_phi_profile), [](auto const& tuple) { return std::get<2>(tuple); });

			// convert vector values to matrices
			Eigen::ArrayX3d pt_profiles(nb_target_samples, 3);
			pt_profiles.col(0) = Eigen::Map<Eigen::ArrayXd>(pt_tau_profile.data(), nb_target_samples);
			pt_profiles.col(1) = Eigen::Map<Eigen::ArrayXd>(pt_kappa_profile.data(), nb_target_samples);
			pt_profiles.col(2) = Eigen::Map<Eigen::ArrayXd>(pt_phi_profile.data(), nb_target_samples);

			// find, for each target point its closest point in source points and its related estimated scale
			std::vector<std::tuple<Point, Scalar, bool>> pairs_source_and_lag{ std::make_tuple(Point(), 0.0, false),  std::make_tuple(Point(), 0.0, false) , std::make_tuple(Point(), 0.0, false) };
			std::multiset<Scalar> max_corr{ 0.0, 0.0, 0.0 };
			//clock_t tstart, tend;
			//tstart = clock();
//#pragma omp parallel for
			for (int v = 0; v < source_gls_profiles.size(); ++v)
			{
				//std::cout << " p_source = " << std::get<0>(source_gls_profiles[v]).pos().transpose() << std::endl;
				const std::vector<std::tuple<Scalar, Scalar, Scalar>> & ps_source_profile = std::get<1>(source_gls_profiles[v]);
				std::vector<Scalar> ps_tau_profile, ps_kappa_profile, ps_phi_profile;

				// get all gls parameter profiles in separate vectors
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_tau_profile), [](auto const& tuple) { return std::get<0>(tuple); });
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_kappa_profile), [](auto const& tuple) { return std::get<1>(tuple); });
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_phi_profile), [](auto const& tuple) { return std::get<2>(tuple); });

				// convert vector values to matrices
				Eigen::ArrayX3d ps_profiles(nb_source_samples, 3);
				ps_profiles.col(0) = Eigen::Map<Eigen::ArrayXd>(ps_tau_profile.data(), nb_source_samples);
				ps_profiles.col(1) = Eigen::Map<Eigen::ArrayXd>(ps_kappa_profile.data(), nb_source_samples);
				ps_profiles.col(2) = Eigen::Map<Eigen::ArrayXd>(ps_phi_profile.data(), nb_source_samples);

				// find relative scale 
					// compute shifts
				std::pair<Scalar, Scalar> pair_shift_Dsigma = compute_optimal_shift(ps_profiles, pt_profiles, W, alpha, ratio);

				//std::cout << "source index : " << v << " , shift : " << pair_shift_Dsigma.first << " , max corr : " << pair_shift_Dsigma.second << std::endl;

					// compare max_Dsigma to max_corr to know if ps and pt belong to one of the 3 matching pairs
				update_maxCorr_and_lag(pair_shift_Dsigma.second, max_corr, pair_shift_Dsigma.first, std::get<0>(source_gls_profiles[v]), pairs_source_and_lag, ratio_corr);
			}

			std::string output_filename = "C:\\Registration\\test_gls_algo\\matching_pairs\\source_matching_pairs_unsorted.txt";
			write_closest_matching_points(std::get<0>(target_gls_profiles[u]), pairs_source_and_lag, output_filename, false);

			// store matching target-source matching pairs with their cost
			for (int k = 0; k < 3; k++)
			{
				const Point & pt = std::get<0>(target_gls_profiles[u]);
				if (std::get<2>(pairs_source_and_lag[k]))
				{
					Scalar cost_pair = target_priority.find(pt)->second * source_priority.find(std::get<0>(pairs_source_and_lag[k]))->second;
					target_source_matchings.push_back(std::make_tuple(pt, std::get<0>(pairs_source_and_lag[k]), std::get<1>(pairs_source_and_lag[k]), cost_pair));
				}
			}

			//tend = clock()-tstart;
			//cout << "time for computing 3 matching pair : " << (float)tend/CLOCKS_PER_SEC << " second(s)." << endl;
		}

		return target_source_matchings;
	}

	else // otherwise match source points on target points
	{
		return compute_3_closest_pairs(target_gls_profiles, source_gls_profiles, ratio, nb_target_samples, nb_source_samples, w, alpha);
	}
}


void PointMap::set_point_profiles_cost(std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>& point_gls_profiles, int nb_samples, int alpha)
{
	size_t nb_points = point_gls_profiles.size();


	for (int k = 0; k<nb_points; k++)
	{
		// compute cost (ie priority)
		Scalar priority = compute_point_priority(std::get<2>(point_gls_profiles[k]), nb_samples, alpha);

		// set profiles
		const std::vector<std::tuple<Scalar, Scalar, Scalar>> & pt_target_profile = std::get<1>(point_gls_profiles[k]);
		std::vector<Scalar> pt_tau_profile, pt_kappa_profile, pt_phi_profile;

			// get all gls parameter profiles in separate vectors
		std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_tau_profile), [](auto const& tuple) { return std::get<0>(tuple); });
		std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_kappa_profile), [](auto const& tuple) { return std::get<1>(tuple); });
		std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_phi_profile), [](auto const& tuple) { return std::get<2>(tuple); });

			// convert vector values to matrices
		Eigen::ArrayX3d pt_profiles(nb_samples, 3);
		pt_profiles.col(0) = Eigen::Map<Eigen::ArrayXd>(pt_tau_profile.data(), nb_samples);
		pt_profiles.col(1) = Eigen::Map<Eigen::ArrayXd>(pt_kappa_profile.data(), nb_samples);
		pt_profiles.col(2) = Eigen::Map<Eigen::ArrayXd>(pt_phi_profile.data(), nb_samples);

		//store values in map
		map_index_point_profiles_cost_.insert(std::make_pair(k, std::make_tuple(std::get<0>(point_gls_profiles[k]), pt_profiles, priority)));
	}
}

std::pair<std::map<int, std::tuple<Point, Eigen::ArrayX3d, Scalar>>::iterator, std::map<int, std::tuple<Point, Eigen::ArrayX3d, Scalar>>::iterator> PointMap::get_iterator_range()
{
	return std::make_pair(map_index_point_profiles_cost_.begin(), map_index_point_profiles_cost_.end());
}

std::map<int, std::tuple<Point, Eigen::ArrayX3d, Scalar>>::iterator PointMap::find_index(int index)
{
	return map_index_point_profiles_cost_.find(index);
}

std::vector<std::tuple<std::pair<int, int>, Scalar, Scalar, Scalar>> compute_symmetric_pairs(PointMap& map_index_profiles_cost_source, PointMap& map_index_profiles_cost_target, Scalar ratio, int k, VectorType w, Scalar alpha, bool cross_check)
{
	Scalar ratio_corr = 0.9;
	Scalar shift_err = 5.0;

	// get W vector   -- slower than [w(0), w(1), w(2)].transpose() ?
	Eigen::Array<Scalar, 1, 3> W;
	W << w[0], w[1], w[2];

	std::map<std::pair<int, int>, std::tuple<Scalar, Scalar, bool>> symmetric_matchings;

	// search for target->source matching
	for (auto it_t  = map_index_profiles_cost_target.get_iterator_range().first; it_t != map_index_profiles_cost_target.get_iterator_range().second; it_t++)
	{
		// intialize max_corr and  pairs_source_and_lag (for k nearest source points)
		std::multiset<Scalar> max_corr;
		std::vector<std::tuple<int, Point, Scalar, Scalar, bool>> pairs_source_lag_and_corr;
		for (int i = 0; i < k; i++)
		{
			max_corr.insert(0.0);
			pairs_source_lag_and_corr.push_back(std::make_tuple(-1, Point(), 0.0, 0.0, false)); 
		}

		// find, for each target point its closest point in source points and its related estimated shift
		for (auto it_s = map_index_profiles_cost_source.get_iterator_range().first; it_s != map_index_profiles_cost_source.get_iterator_range().second; it_s++)
		{
			// compute shifts
			std::pair<Scalar, Scalar> pair_shift_Dsigma = compute_optimal_shift(std::get<1>(it_s->second), std::get<1>(it_t->second), W, alpha, ratio);

			// compare max_Dsigma to max_corr to know if source belong to one of the k matching pairs
			update_maxCorr_and_lag(pair_shift_Dsigma.second, max_corr, pair_shift_Dsigma.first, it_s-> first, std::get<0>(it_s->second), pairs_source_lag_and_corr, ratio_corr);

		}

		std::string output_filename = "C:\\Registration\\test_gls_algo\\matching_pairs\\source_matching_pairs_unsorted.txt";
		write_closest_matching_points(std::get<0>(it_t->second), pairs_source_lag_and_corr, output_filename, false);

		// store matching target-source matching pairs with their cost
		for (int i = 0; i < k; i++)
		{
			if (std::get<4>(pairs_source_lag_and_corr[i]))
				symmetric_matchings.insert(std::make_pair(std::make_pair(it_t->first, std::get<0>(pairs_source_lag_and_corr[i])), std::make_tuple(std::get<2>(pairs_source_lag_and_corr[i]), std::get<3>(pairs_source_lag_and_corr[i]), false)));
				
		}

	}

	if (cross_check)
	{
		std::map<int, std::tuple<Point, Eigen::ArrayX3d, Scalar>> symmetric_source_pairs;
		// get source points needed for cross_check
		for (auto it = symmetric_matchings.begin(); it!= symmetric_matchings.end(); it++)
				symmetric_source_pairs.insert(*(map_index_profiles_cost_source.find_index(it->first.second)));
			

		// search for symmetric source->target matching
		for (auto it_s = symmetric_source_pairs.begin(); it_s != symmetric_source_pairs.end(); it_s++)
		{
			// intialize max_corr and  pairs_target_and_lag (for k nearest source points)
			std::multiset<Scalar> max_corr;
			std::vector<std::tuple<int, Point, Scalar, Scalar, bool>> pairs_target_lag_and_corr;
			for (int i = 0; i < k; i++)
			{
				max_corr.insert(0.0);
				pairs_target_lag_and_corr.push_back(std::make_tuple(-1, Point(), 0.0, 0.0, false));
			}

			// find, for each symmetric source point its closest point in target and its related estimated shift
			for (auto it_t = map_index_profiles_cost_target.get_iterator_range().first; it_t != map_index_profiles_cost_target.get_iterator_range().second; it_t++)
			{
				// compute shifts
				std::pair<Scalar, Scalar> pair_shift_Dsigma = compute_optimal_shift(std::get<1>(it_t->second), std::get<1>(it_s->second), W, alpha, ratio);

				// compare max_Dsigma to max_corr to know if target point belongs to one of the k matching pairs
				update_maxCorr_and_lag(pair_shift_Dsigma.second, max_corr, pair_shift_Dsigma.first, it_t->first, std::get<0>(it_t->second), pairs_target_lag_and_corr, ratio_corr);

			}

			// check if source->target matchings are compatible with target->source matchings
			for (int i = 0; i < k; i++)
			{
				if (std::get<4>(pairs_target_lag_and_corr[i]))
				{
					auto it = symmetric_matchings.find(std::make_pair(std::get<0>(pairs_target_lag_and_corr[i]), it_s->first));
					if (it != symmetric_matchings.end())
					{
						// check if scales values are close
						Scalar scale_s = std::get<2>(pairs_target_lag_and_corr[i]);
						Scalar scale_t = std::get<0>(it->second);
						if (abs(scale_s + scale_t) <= shift_err)
							(it->second) = std::make_tuple(std::get<0>(it->second), std::get<1>(it->second), true);	 // indicate that a symmetric matching pair has been found
					}

				}
			}

		}


		// collect matching pairs
		int pos = 0;
		auto it_begin = symmetric_matchings.begin();
		auto it = symmetric_matchings.begin();
		while (it != symmetric_matchings.end())
		{
			if (!std::get<2>(it->second))
				symmetric_matchings.erase(it);
			else
				pos += 1;

			it = std::next(it_begin, pos);
		}

	}

	std::vector<std::tuple<std::pair<int, int>, Scalar, Scalar, Scalar>> target_source_matchings;
	for (auto it = symmetric_matchings.begin(); it != symmetric_matchings.end(); it++)
	{
		Scalar lag = std::get<0>(it->second);
		Scalar cost = std::get<2>((map_index_profiles_cost_target.find_index(it->first.first))->second)*std::get<2>((map_index_profiles_cost_source.find_index(it->first.second))->second);
		Scalar corr = std::get<1>(it->second);
		target_source_matchings.push_back(std::make_tuple(it->first, lag, cost, corr));

	}

	return target_source_matchings;
		
}


std::pair<Scalar, Scalar> compute_optimal_shift(Eigen::ArrayX3d& source_profiles, Eigen::ArrayX3d& target_profiles, Eigen::Array<Scalar, 1, 3>& W, Scalar alpha, Scalar ratio, bool reverse)
{
	int p = source_profiles.rows();
	int q = target_profiles.rows();

	Scalar opt_Dsigma = 0.0;
	Scalar opt_shift = 0.0;

	if (p >= q)
	{
		int min_shift = -floor((1.0 - ratio)*q);
		int max_shift = (p - q) + floor((1.0 - ratio)*(Scalar)q);

		for (int shift = min_shift; shift < max_shift + 1; shift ++)
		{
			Scalar Dsigma = 0.0;
			if (shift >= 0 && shift <= p - q)  // regular case
			{
				Eigen::ArrayX3d diff(q, 3);
				diff.block(0, 0, q, 3) = source_profiles.block(shift, 0, q, 3) - target_profiles;
				//diff.block(0, 0, q, 3) = (1.0/ source_profiles.block(shift, 0, q, 3).abs().colwise().maxCoeff())*source_profiles.block(shift, 0, q, 3) - (1.0/ target_profiles.abs().colwise().maxCoeff())*target_profiles;
				Eigen::Array<Scalar, 1, 3> den_W = get_max(source_profiles.block(shift, 0, q, 3).abs().colwise().maxCoeff(), target_profiles.abs().colwise().maxCoeff());
				//den_W(2) = 1.0;
				W = 1.0 / den_W;
				//W = 1.0/(diff.abs().colwise().maxCoeff());
				//std::cout << "W : " << W << std::endl;
				Eigen::Array<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> diff2 = (diff*diff).rowwise()*W;
				Eigen::Array<Scalar, Eigen::Dynamic, 1> diss_over_scales = diff2.rowwise().sum();
				Eigen::ArrayXd sigma = 1.0 - (alpha*diss_over_scales).tanh();
				Dsigma = 1.0 / (Scalar)(q) * sigma.sum();
				//std::cout << "Dsigma : " << Dsigma << std::endl;
			}
			else  // extremity case
			{
				if (shift < 0)  // left extremity
				{
					Eigen::ArrayX3d diff(q - abs(shift), 3);
					diff.block(0, 0, q - abs(shift), 3) = source_profiles.block(0, 0, q - abs(shift), 3) - target_profiles.block(abs(shift), 0, q - abs(shift), 3);
					//diff.block(0, 0, q - abs(shift), 3) = (1.0/ source_profiles.block(0, 0, q - abs(shift), 3).abs().colwise().maxCoeff())*source_profiles.block(0, 0, q - abs(shift), 3) - (1.0/ target_profiles.block(abs(shift), 0, q - abs(shift), 3).abs().colwise().maxCoeff())*target_profiles.block(abs(shift), 0, q - abs(shift), 3);
					Eigen::Array<Scalar, 1, 3> den_W = get_max(source_profiles.block(0, 0, q - abs(shift), 3).abs().colwise().maxCoeff(), target_profiles.block(abs(shift), 0, q - abs(shift), 3).abs().colwise().maxCoeff());
					//den_W(2) = 1.0;
					W = 1.0 / den_W;
					//W = 1.0 / (diff.abs().colwise().maxCoeff());
					//std::cout << "W : " << W << std::endl;
					Eigen::Array<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> diff2 = (diff*diff).rowwise()*W;
					Eigen::Array<Scalar, Eigen::Dynamic, 1> diss_over_scales = (diff2).rowwise().sum();
					Eigen::ArrayXd sigma = 1.0 - (alpha*diss_over_scales).tanh();
					Dsigma = 1.0 / (Scalar)(q - abs(shift)) * sigma.sum();
					//std::cout << "Dsigma : " << Dsigma << std::endl;
				}

				else  // right extremity
				{
					Eigen::ArrayX3d diff(p - shift, 3);
					diff.block(0, 0, p - shift, 3) = source_profiles.block(shift, 0, p - shift, 3) - target_profiles.block(0, 0, p - shift, 3);
					//diff.block(0, 0, q - abs(shift), 3) = (1.0/ source_profiles.block(shift, 0, p - shift, 3).abs().colwise().maxCoeff())*source_profiles.block(shift, 0, p - shift, 3) - (1.0 / target_profiles.block(0, 0, p - shift, 3).abs().colwise().maxCoeff())*target_profiles.block(0, 0, p - shift, 3);
					Eigen::Array<Scalar, 1, 3> den_W = get_max(source_profiles.block(shift, 0, p - shift, 3).abs().colwise().maxCoeff(), target_profiles.block(0, 0, p - shift, 3).abs().colwise().maxCoeff());
					//den_W(2) = 1.0;
					W = 1.0 / den_W;
					//W = 1.0/ (diff.abs().colwise().maxCoeff());
					//std::cout << "W : " << W << std::endl;
					Eigen::Array<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> diff2 = (diff*diff).rowwise()*W;
					Eigen::Array<Scalar, Eigen::Dynamic, 1> diss_over_scales = (diff2).rowwise().sum();
					Eigen::ArrayXd sigma = 1.0 - (alpha*diss_over_scales).tanh();
					Dsigma = 1.0 / (Scalar)(p - shift) * sigma.sum();
					//std::cout << "Dsigma : " << Dsigma << std::endl;
				}
			}

			if (Dsigma > opt_Dsigma)
			{
				opt_Dsigma = Dsigma;
				opt_shift = shift;
			}
		}

		if (reverse)
			opt_shift = -opt_shift;

		return std::make_pair(opt_shift, opt_Dsigma);

	}

	else
	{
		return compute_optimal_shift(target_profiles, source_profiles, W, alpha, ratio, true);
	}

}



std::pair<Scalar, std::vector<std::tuple<Point, Point, Scalar>>> relative_scale_estimation(std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& source_gls_profiles, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& target_gls_profiles, Scalar ratio, int nb_samples, VectorType w, Scalar alpha)
{
	size_t nb_source_points = source_gls_profiles.size();
	size_t nb_target_points = target_gls_profiles.size();

	if (nb_source_points >= nb_target_points)   // if nb_source_points >= nb_target_points, then match target points on source points
	{
		std::vector<std::tuple<Point, Point, Scalar>> target_source_scale_matchings;

		for (std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>::iterator it_t = target_gls_profiles.begin(); it_t != target_gls_profiles.end(); ++it_t)
		{
			std::vector<std::tuple<Scalar, Scalar, Scalar>> pt_target_profile = it_t->second;
			std::vector<Scalar> pt_tau_profile, pt_kappa_profile, pt_phi_profile;

			// get all gls parameter profiles in separate vectors
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_tau_profile), [](auto const& tuple) { return std::get<0>(tuple); });
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_kappa_profile), [](auto const& tuple) { return std::get<1>(tuple); });
			std::transform(std::begin(pt_target_profile), std::end(pt_target_profile), std::back_inserter(pt_phi_profile), [](auto const& tuple) { return std::get<2>(tuple); });

			// convert vector values to matrices
			Eigen::ArrayX3d pt_profiles = Eigen::ArrayX3d::Zero(nb_samples, 3);
			pt_profiles.col(0) = Eigen::Map<Eigen::ArrayXd>(pt_tau_profile.data(), nb_samples);
			pt_profiles.col(1) = Eigen::Map<Eigen::ArrayXd>(pt_kappa_profile.data(), nb_samples);
			pt_profiles.col(2) = Eigen::Map<Eigen::ArrayXd>(pt_phi_profile.data(), nb_samples);

			/*std::cout << "pt profile : " << std::endl;
			std::cout << pt_profiles << std::endl;*/

			// find, for each target point its closest point in source points and its related estimated scale
			std::pair<Point, Scalar> source_pair_and_scale;
			Scalar max_corr =  0.0;

			for (std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>::iterator it_s = source_gls_profiles.begin(); it_s != source_gls_profiles.end(); ++it_s)
			{
				std::vector<std::tuple<Scalar, Scalar, Scalar>> ps_source_profile = it_s->second;
				std::vector<Scalar> ps_tau_profile, ps_kappa_profile, ps_phi_profile;

				// get all gls parameter profiles in separate vectors
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_tau_profile), [](auto const& tuple) { return std::get<0>(tuple); });
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_kappa_profile), [](auto const& tuple) { return std::get<1>(tuple); });
				std::transform(std::begin(ps_source_profile), std::end(ps_source_profile), std::back_inserter(ps_phi_profile), [](auto const& tuple) { return std::get<2>(tuple); });

				// convert vector values to matrices
				Eigen::ArrayX3d ps_profiles = Eigen::ArrayX3d::Zero(nb_samples, 3);
				ps_profiles.col(0) = Eigen::Map<Eigen::ArrayXd>(ps_tau_profile.data(), nb_samples);
				ps_profiles.col(1) = Eigen::Map<Eigen::ArrayXd>(ps_kappa_profile.data(), nb_samples);
				ps_profiles.col(2) = Eigen::Map<Eigen::ArrayXd>(ps_phi_profile.data(), nb_samples);

				/*std::cout << "ps profile : " << std::endl;
				std::cout << ps_profiles << std::endl;*/


				// find relative scale 
				Scalar max_Dsigma = 0.0;
				Scalar lag = 0.0;
					// compute shifts
				Scalar max_shift = nb_samples - floor(nb_samples*ratio);
						// right shifts
				for (int shift = 0; shift < max_shift -2; shift++)
				{	 
					Eigen::ArrayX3d diff = Eigen::ArrayX3d::Zero(nb_samples-shift, 3);

					diff.block(0, 0, nb_samples - shift, 3) = ps_profiles.block(0, 0, nb_samples-shift, 3) - pt_profiles.block(shift, 0, nb_samples, 3);
					/*std::cout << "diff size : " << diff.rows() << std::endl;
					std::cout << "diff : " << std::endl;
					std::cout << diff << std::endl;*/

					Eigen::ArrayX3d diff2 = (diff*diff).rowwise()*Eigen::Array3d(w[0], w[1], w[2]).transpose();
					/*std::cout << "diff2 size : " << diff2.rows() << std::endl;
					std::cout << "diff2 : " << std::endl;
					std::cout << diff2 << std::endl;*/

					Eigen::ArrayXd diss_over_scales = (diff2).rowwise().sum();
					/*std::cout << "diss : " << std::endl;
					std::cout << diss_over_scales << std::endl;*/

					Eigen::ArrayXd sigma = 1.0 - (alpha*diss_over_scales).tanh();
					/*std::cout << "sigma : " << std::endl;
					std::cout << sigma << std::endl;
					std::cout << "sigma size : " << sigma.size() << std::endl;*/

					Scalar Dsigma = 1.0/(Scalar)(nb_samples-1-shift) * sigma.sum();
					/*std::cout << "Dsigma : " << std::endl;
					std::cout << Dsigma << std::endl;*/

					if (Dsigma > max_Dsigma)
					{
						max_Dsigma = Dsigma;
						lag = -shift;
					}
				}
						// left shifts
				for (int shift = 0; shift < max_shift - 2; shift++)
				{
					Eigen::ArrayX3d diff = Eigen::ArrayX3d::Zero(nb_samples-shift, 3);

					diff.block(0, 0, nb_samples - shift, 3) = ps_profiles.block(shift, 0, nb_samples, 3) - pt_profiles.block(0, 0, nb_samples - shift, 3);
					Eigen::ArrayX3d diff2 = (diff*diff).rowwise()*Eigen::Array3d(w[0], w[1], w[2]).transpose();
					Eigen::ArrayXd diss_over_scales = (diff2).rowwise().sum();
					Eigen::ArrayXd sigma = 1.0 - (alpha*diss_over_scales).tanh();
					Scalar Dsigma = 1.0 / (Scalar)(nb_samples - 1 - shift) * sigma.sum();

					if (Dsigma > max_Dsigma)
					{
						max_Dsigma = Dsigma;
						lag = shift;
					}
				}

				// compare max_Dsigma to max_corr to know if ps and pt are matching pair
				if (max_Dsigma > max_corr)
				{
					max_corr = max_Dsigma;
					source_pair_and_scale = std::make_pair(it_s->first, lag);
				}

			}

			// store matching target-source matching pairs
			target_source_scale_matchings.push_back(std::make_tuple(it_t->first, source_pair_and_scale.first, source_pair_and_scale.second));	
		}

		// compute relative scale
		std::vector<Scalar> relative_scales;
		std::transform(std::begin(target_source_scale_matchings), std::end(target_source_scale_matchings), std::back_inserter(relative_scales), [](auto const& tuple) { return std::get<2>(tuple); });
		//return std::make_pair(std::accumulate(relative_scales.begin(), relative_scales.end(), 0.0) / relative_scales.size(), target_source_scale_matchings);  // if we take relative scale as the average of best scales for target->source matching pairs
		return std::make_pair(Median(relative_scales.begin(), relative_scales.end()), target_source_scale_matchings);  // if we take relative scale as the median of best scales for target->source matching pairs
	}

	else // otherwise match source points on target points
	{
		return relative_scale_estimation(target_gls_profiles, source_gls_profiles, ratio, nb_samples, w, alpha);
	}
}
