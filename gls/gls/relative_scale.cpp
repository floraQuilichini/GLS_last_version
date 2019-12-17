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
#include "matching_prioritization.h"

using namespace std;


Scalar Median(std::vector<Scalar>::iterator begin, std::vector<Scalar>::iterator end)
{
	const auto size = std::distance(begin, end);
	std::nth_element(begin, begin + size / 2, end);
	return *std::next(begin, size / 2);
}


void update_maxCorr_and_lag(Scalar max_Dsigma, std::multiset<Scalar>& max_corr, Scalar lag, Point point, std::vector<std::pair<Point, Scalar>>& pairs_point_and_scale)
{
	std::multiset<Scalar>::iterator it = max_corr.insert(max_Dsigma);
	int dist = std::distance(max_corr.begin(), it);
	max_corr.erase(max_corr.begin());

	if (dist > 0)
		pairs_point_and_scale[dist-1] = std::make_pair(point, lag);
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
			std::vector<std::pair<Point, Scalar>> pairs_source_and_scale{ std::make_pair(Point(), 0.0),  std::make_pair(Point(), 0.0) , std::make_pair(Point(), 0.0) };
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
				update_maxCorr_and_lag(max_Dsigma, max_corr, lag, it_s->first, pairs_source_and_scale);

			}

			// store matching target-source matching pairs
			for (int k =0; k<3; k++ )
				target_source_scale_matchings.push_back(std::make_tuple(it_t->first, pairs_source_and_scale[k].first, pairs_source_and_scale[k].second));
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

	if (nb_source_points >= nb_target_points)   // if nb_source_points >= nb_target_points, then match target points on source points
	{
		std::vector<std::tuple<Point, Point, Scalar, Scalar>> target_source_matchings;

		//compute point priority
		std::map<Point, Scalar> target_priority, source_priority;
		for (int k=0; k<nb_target_points; k++)
		{
			Scalar priority = compute_point_priority(std::get<2>(target_gls_profiles[k]), nb_samples, alpha);
			target_priority.insert(std::make_pair(std::get<0>(target_gls_profiles[k]), priority));
		}

		for (int k = 0; k<nb_source_points; k++)
		{
			Scalar priority = compute_point_priority(std::get<2>(source_gls_profiles[k]), nb_samples, alpha);
			source_priority.insert(std::make_pair(std::get<0>(source_gls_profiles[k]), priority));
		}

		// search for target->source matching
#pragma omp parallel for
		for (int u = 0 ; u < target_gls_profiles.size() ; ++u)
		// for (std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>::iterator it_t = target_gls_profiles.begin(); it_t != target_gls_profiles.end(); ++it_t)
		{
			const std::vector<std::tuple<Scalar, Scalar, Scalar>> & pt_target_profile = std::get<1>(target_gls_profiles[u]);
			// std::vector<std::tuple<Scalar, Scalar, Scalar>> pt_target_profile = std::get<1>(*it_t);
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
			std::vector<std::pair<Point, Scalar>> pairs_source_and_scale{ std::make_pair(Point(), 0.0),  std::make_pair(Point(), 0.0) , std::make_pair(Point(), 0.0) };
			std::multiset<Scalar> max_corr{ 0.0, 0.0, 0.0 };
			clock_t tstart, tend;
			tstart = clock();
			for (std::vector<std::tuple<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>, std::vector<Scalar>>>::iterator it_s = source_gls_profiles.begin(); it_s != source_gls_profiles.end(); ++it_s)
			{
				const std::vector<std::tuple<Scalar, Scalar, Scalar>> & ps_source_profile = std::get<1>(*it_s);
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
				Scalar max_Dsigma = 0.0;
				Scalar lag = 0.0;
				// compute shifts
				Scalar max_shift = nb_samples - floor(nb_samples*ratio);
				// right shifts
				for (int shift = 0; shift < max_shift - 2; shift++)
				{
					Eigen::ArrayX3d diff = Eigen::ArrayX3d::Zero(nb_samples - shift, 3);
					//std::cout << "diff nb rows : " << diff.rows() << " diff nb cols : " << diff.cols() << std::endl;
					diff.block(0, 0, nb_samples - shift, 3) = ps_profiles.block(0, 0, nb_samples - shift, 3) - pt_profiles.block(shift, 0, nb_samples - shift, 3);
					//std::cout << " diff size : " << diff.size() << std::endl;
					Eigen::ArrayX3d diff2 = (diff*diff).rowwise()*Eigen::Array3d(w[0], w[1], w[2]).transpose();
					//std::cout << " diff2 : " << diff2.size() << std::endl;
					Eigen::ArrayXd diss_over_scales = (diff2).rowwise().sum();
					//std::cout << " diss_over_scales : " << diss_over_scales.size() << std::endl;
					Eigen::ArrayXd sigma = 1.0 - (alpha*diss_over_scales).tanh();
					//std::cout << " sigma : " << sigma.size() << std::endl;
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

					diff.block(0, 0, nb_samples - shift, 3) = ps_profiles.block(shift, 0, nb_samples - shift, 3) - pt_profiles.block(0, 0, nb_samples - shift, 3);
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
				update_maxCorr_and_lag(max_Dsigma, max_corr, lag, std::get<0>(*it_s), pairs_source_and_scale);

			}

			// store matching target-source matching pairs with their cost
			for (int k = 0; k < 3; k++)
			{
				const Point & pt = std::get<0>(target_gls_profiles[u]);
				Scalar cost_pair = target_priority.find(pt)->second * source_priority.find(pairs_source_and_scale[k].first)->second;
				target_source_matchings.push_back(std::make_tuple(pt, pairs_source_and_scale[k].first, pairs_source_and_scale[k].second, cost_pair));
				//Scalar cost_pair = target_priority.find(std::get<0>(*it_t))->second * source_priority.find(pairs_source_and_scale[k].first)->second;
				//target_source_matchings.push_back(std::make_tuple(std::get<0>(*it_t), pairs_source_and_scale[k].first, pairs_source_and_scale[k].second, cost_pair));
			}

			tend = clock()-tstart;
			cout << "time for computing 3 matching pair : " << (float)tend/CLOCKS_PER_SEC << " second(s)." << endl;
		}

		return target_source_matchings;
	}

	else // otherwise match source points on target points
	{
		return compute_3_closest_pairs(target_gls_profiles, source_gls_profiles, ratio, nb_samples, w, alpha);
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
