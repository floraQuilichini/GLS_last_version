#include "stdafx.h"

#include <cmath>
#include <algorithm>
#include <vector>
#include <tuple>
#include <string>
#include <iterator>
#include <chrono>
#include <random>
#include "Patate/grenaille.h"
#include "Eigen/Eigen"
#include "IO.h"
#include "ransac_scheme.h"


Scalar RansacScheme::compute_points_dist(Point& p1, Point& p2)
{
	return sqrt(((p1.pos().transpose() - p2.pos().transpose())*(p1.pos() - p2.pos())).sum());
}


std::pair<Scalar, int> RansacScheme::compute_min(std::vector<std::pair<Scalar, int>>& vec, int index_to_remove)
{
	std::pair<Scalar, int> min_index;
	if (index_to_remove > -1)
	{

		if (vec[0].second == index_to_remove)
			min_index = std::make_pair(vec[1].first, 0);
		else
			min_index = std::make_pair(vec[0].first, 0);
		for (int i = 0; i < (int)vec.size(); i++)
		{
			if (vec[i].first < min_index.first && vec[i].second != index_to_remove)
				min_index = std::make_pair(vec[i].first, i);
		}
	}

	else
	{
		min_index = std::make_pair(vec[0].first, 0);
		for (int i = 0; i < (int)vec.size(); i++)
		{
			if (vec[i].first < min_index.first)
				min_index = std::make_pair(vec[i].first, i);
		}
	}

	return min_index;
}


void RansacScheme::compute_target_optimal_costs()
{
	auto it_range = target_pointMap_->get_iterator_range();
	max_cost_ = std::pow(std::get<2>((it_range.first)->second), 2.0);
	min_cost_ = max_cost_;
	for (auto it = it_range.first; it != it_range.second; it++)
	{
		Scalar cost = std::pow(std::get<2>(it->second), 2.0);
		if (cost > max_cost_)
			max_cost_ = cost;
		else
		{
			if (cost < min_cost_)
				min_cost_ = cost;
		}
	}
}



/*void RansacScheme::find_2_farthest_pairs(RansacScheme::triplet& triplet)
{
	Scalar dist2 = 0.0;
	Scalar dist3 = 0.0;
	std::tuple<Point, Point, Scalar, Scalar> pair2, pair3;
	Point p1 = std::get<0>(triplet.pair1);
	queue_copy_ = *queue_ptr_;

	while(!queue_copy_.empty())
	{
		std::tuple<Point, Point, Scalar, Scalar> pair = queue_copy_.top();
		queue_copy_.pop();
		if (is_q_unique(triplet, pair))
		{
			Point p = std::get<0>(pair);
			Scalar dist = compute_points_dist(p1, p);
			if (dist > dist2)
			{
				triplet.pair3 = triplet.pair2;
				dist3 = dist2;
				triplet.pair2 = pair;
				dist2 = dist;
			}
			else
			{
				if (dist > dist3)
				{
					triplet.pair3 = pair;
					dist3 = dist;
				}
			}
		}
	}

}*/

void RansacScheme::find_2_farthest_pairs(RansacScheme::triplet& triplet)
{
	// get second pair
	triplet.pair2 = find_farthest_pair(triplet.pair1);

	// get third pair
	triplet.pair3 = find_pair_to_maximize_triangle_surface(triplet.pair1, triplet.pair2);

}



std::tuple<int, int, Scalar, Scalar> RansacScheme::find_farthest_pair(std::tuple<int, int, Scalar, Scalar>& pair)
{
	Scalar dist_ref = 0.0;
	std::tuple<int, int, Scalar, Scalar> farthest_pair;
	Point p = std::get<0>(target_pointMap_->find_index(std::get<0>(pair))->second);
	queue_copy_ = *queue_ptr_;

	while (!queue_copy_.empty())
	{
		std::tuple<int, int, Scalar, Scalar> query_pair = queue_copy_.top();
		queue_copy_.pop();
		Point p_query = std::get<0>(target_pointMap_->find_index(std::get<0>(query_pair))->second);
		Scalar dist = compute_points_dist(p, p_query);
		if (dist > dist_ref && !(std::get<0>(source_pointMap_->find_index(std::get<1>(query_pair))->second).pos() == std::get<0>(source_pointMap_->find_index(std::get<1>(pair))->second).pos()))
		{
			dist_ref = dist_ref;
			farthest_pair = query_pair;
		}
	}

	return farthest_pair;
}

std::tuple<int, int, Scalar, Scalar> RansacScheme::find_pair_to_maximize_triangle_surface(std::tuple<int, int, Scalar, Scalar>& pair1, std::tuple<int, int, Scalar, Scalar>& pair2)
{
	Point p1 = std::get<0>(target_pointMap_->find_index(std::get<0>(pair1))->second);
	Point p2 = std::get<0>(target_pointMap_->find_index(std::get<0>(pair2))->second);
	std::tuple<int, int, Scalar, Scalar> pair3;
	Scalar max_area = 0.0;
	Scalar a = compute_points_dist(p1, p2);
	queue_copy_ = *queue_ptr_;

	while (!queue_copy_.empty())
	{
		std::tuple<int, int, Scalar, Scalar> pair = queue_copy_.top();
		queue_copy_.pop();
		if (!(std::get<0>(source_pointMap_->find_index(std::get<1>(pair1))->second).pos() == std::get<0>(source_pointMap_->find_index(std::get<1>(pair))->second).pos() || std::get<0>(source_pointMap_->find_index(std::get<1>(pair2))->second).pos() == std::get<0>(source_pointMap_->find_index(std::get<1>(pair))->second).pos()))
		{
			// compute area of the triangle formed by the point std::get<1>(pair1), std::get<1>(pair2), std::get<1>(pair)
			Scalar b = compute_points_dist(p1, std::get<0>(target_pointMap_->find_index(std::get<0>(pair))->second));
			Scalar c = compute_points_dist(p2, std::get<0>(target_pointMap_->find_index(std::get<0>(pair))->second));
			Scalar area = sqrt((a+b+c)*(b+c-a)*(a+c-b)*(a+b-c))/4.0;
			if (area > max_area)
			{
				max_area = area;
				pair3 = pair;
			}
		}

	}

	return pair3;

}


void RansacScheme::find_2_farthest_pairs_with_respect_to_geomVar(RansacScheme::triplet& triplet, Scalar bbox_diag, Scalar lambda)
{
	queue_copy_ = *queue_ptr_;

	// get second pair
	Scalar max_metric_value = 0.0;
	Point& pt1 = std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair1))->second);
	Scalar priority_pt1 = std::get<2>(target_pointMap_->find_index(std::get<0>(triplet.pair1))->second);
	Point& ps1 = std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair1))->second);
	while (!queue_copy_.empty())
	{
		std::tuple<int, int, Scalar, Scalar> pair = queue_copy_.top();
		queue_copy_.pop();
		Point& query_s_point = std::get<0>(source_pointMap_->find_index(std::get<1>(pair))->second);
		Point& query_t_point = std::get<0>(target_pointMap_->find_index(std::get<0>(pair))->second);
		Scalar priority_query_t_point = std::get<2>(target_pointMap_->find_index(std::get<0>(pair))->second);
		Scalar distL2 = compute_points_dist(pt1, query_t_point);
		Scalar cost = priority_pt1* priority_query_t_point;
		Scalar metric_value = distL2/bbox_diag + lambda*((max_cost_ - cost)/(max_cost_ - min_cost_));
		if (!(ps1.pos() == query_s_point.pos() || pt1.pos() == query_t_point.pos())) // if pairs are different
		{
			if (metric_value > max_metric_value)
			{
				triplet.pair2 = pair;
				max_metric_value = metric_value;
			}
		}
	}

	// get third pair
	queue_copy_ = *queue_ptr_;
	max_metric_value = 0.0;
	Point& pt2 = std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair2))->second);
	Scalar priority_pt2 = std::get<2>(target_pointMap_->find_index(std::get<0>(triplet.pair2))->second);
	Point& ps2 = std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair2))->second);

	while (!queue_copy_.empty())
	{
		std::tuple<int, int, Scalar, Scalar> pair = queue_copy_.top();
		queue_copy_.pop();
		Point& query_s_point = std::get<0>(source_pointMap_->find_index(std::get<1>(pair))->second);
		Point& query_t_point = std::get<0>(target_pointMap_->find_index(std::get<0>(pair))->second);
		Scalar priority_query_t_point = std::get<2>(target_pointMap_->find_index(std::get<0>(pair))->second);
		Scalar distL2 = 0.5*(compute_points_dist(pt1, query_t_point) + compute_points_dist(pt2, query_t_point));
		Scalar cost = 0.5*(priority_pt1 * priority_query_t_point + priority_pt2 * priority_query_t_point);
		Scalar metric_value = distL2 / bbox_diag + lambda*((max_cost_ - cost)/(max_cost_ - min_cost_));

		if (!(ps1.pos() == query_s_point.pos() || ps2.pos() == query_s_point.pos() || pt1.pos() == query_t_point.pos() || pt2.pos() == query_t_point.pos())) // if pairs are different
		{
			if (metric_value > max_metric_value)
			{
				triplet.pair3 = pair;
				max_metric_value = metric_value;
			}
		}
	}

	// empty queue
	while (!queue_copy_.empty())
		queue_copy_.pop();

}


void RansacScheme::find_2_farthest_pairs_with_respect_to_geomVar(RansacScheme::triplet& triplet, Scalar bbox_diag_ratio)
{
	queue_copy_ = *queue_ptr_;

	// get second pair
	while (!queue_copy_.empty())
	{
		std::tuple<int, int, Scalar, Scalar> pair = queue_copy_.top();
		queue_copy_.pop();
		if (compute_points_dist(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair1))->second), std::get<0>(target_pointMap_->find_index(std::get<0>(pair))->second)) > bbox_diag_ratio && !(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair1))->second).pos() == std::get<0>(source_pointMap_->find_index(std::get<1>(pair))->second).pos()))
		{
			triplet.pair2 = pair;
			break;
		}
	}
		// empty queue
	while (!queue_copy_.empty())
		queue_copy_.pop();


	// get third pair
	queue_copy_ = *queue_ptr_;
	Scalar a = compute_points_dist(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair1))->second), std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair2))->second));
	Scalar max_area = sqrt(3.0) / 4.0* bbox_diag_ratio*bbox_diag_ratio;
	while (!queue_copy_.empty())
	{
		std::tuple<int, int, Scalar, Scalar> pair = queue_copy_.top();
		queue_copy_.pop();
		if (!(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair1))->second).pos() == std::get<0>(source_pointMap_->find_index(std::get<1>(pair))->second).pos() || std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair2))->second).pos() == std::get<0>(source_pointMap_->find_index(std::get<1>(pair))->second).pos()))
		{
			// compute area of the triangle formed by the point std::get<1>(pair1), std::get<1>(pair2), std::get<1>(pair)
			Scalar b = compute_points_dist(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair1))->second), std::get<0>(target_pointMap_->find_index(std::get<0>(pair))->second));
			Scalar c = compute_points_dist(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair2))->second), std::get<0>(target_pointMap_->find_index(std::get<0>(pair))->second));
			Scalar area = sqrt((a + b + c)*(b + c - a)*(a + c - b)*(a + b - c)) / 4.0;
			if (area > max_area)
			{
				triplet.pair3 = pair;
				break;
			}
		}
	}

	// empty queue
	while (!queue_copy_.empty())
		queue_copy_.pop();

}

RansacScheme::RansacScheme(pair_priority_queue& pair_prioritized)
{
	queue_ptr_ = pair_prioritized.get_queue_ptr();
	compute_target_optimal_costs();
	//std::cout << "min cost : " << min_cost_ << std::endl;
	//std::cout << "max cost : " << max_cost_ << std::endl;
}
RansacScheme::RansacScheme(pair_priority_queue& pair_prioritized, PointMap* source_pointMap, PointMap* target_pointMap)
{
	queue_ptr_ = pair_prioritized.get_queue_ptr();
	source_pointMap_ = source_pointMap;
	target_pointMap_ = target_pointMap;
	compute_target_optimal_costs();
	//std::cout << "min cost : " << min_cost_ << std::endl;
	//std::cout << "max cost : " << max_cost_ << std::endl;
}

RansacScheme::triplet RansacScheme::pop_triplet()
{
	RansacScheme::triplet triplet;
	bool is_triplet_complete = false;

	while (queue_ptr_->size() > 2 && !is_triplet_complete)
	{
		// get first pair
		triplet.pair1 = queue_ptr_->top();
		queue_ptr_->pop();

		//get two other pairs
		queue_copy_ = *queue_ptr_;
		int counter = 1;
		while (!queue_copy_.empty())
		{
			std::tuple<int, int, Scalar, Scalar> pair = queue_copy_.top();
			queue_copy_.pop();
			if (!(std::get<0>(target_pointMap_->find_index(std::get<0>(pair))->second) == std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair1))->second) || std::get<0>(source_pointMap_->find_index(std::get<1>(pair))->second) == std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair1))->second)))
			{
				if (counter == 1)
				{
					triplet.pair2 = pair;
					counter += 1;
				}
				if (counter == 2)
				{
					if (!(std::get<0>(target_pointMap_->find_index(std::get<0>(pair))->second) == std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair2))->second) || std::get<0>(source_pointMap_->find_index(std::get<1>(pair))->second) == std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair2))->second)))
					{
						triplet.pair3 = pair;
						is_triplet_complete = true;
						break;
					}
				}
				
			}
		}

		// empty queue copy
		while (!queue_copy_.empty())
			queue_copy_.pop();
	}

	return triplet;
}


RansacScheme::triplet RansacScheme::pick_triplet(const int ind1, const int ind2, const int ind3)
{
	queue_copy_ = *queue_ptr_;
	RansacScheme::triplet triplet;
	int counter = 0;
	while (!queue_copy_.empty())
	{
		if (counter + 1 == ind1)
			triplet.pair1 = queue_copy_.top();
		else
		{
			if (counter + 1 == ind2)
				triplet.pair2 = queue_copy_.top();
			else
			{
				if (counter + 1 == ind3)
					triplet.pair3 = queue_copy_.top();
			}
		}

		queue_copy_.pop();
		counter += 1;
	 }

	return triplet;
}


std::vector<RansacScheme::triplet> RansacScheme::pick_triplets_between_k_pairs(RansacScheme::pqueue& kpairs, int k, int nb_iter)
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_int_distribution<int> distribution(0, k - 1);

	std::vector<triplet> vec_triplets;
	for (int i = 0; i < nb_iter; i++)
	{
		// get random triplet of integers
		std::set<int> indices_in_ascending_order{ k + 1 };
		int counter = 0;
		size_t size_set = 0;
		do
		{
			//int index = rand() % k;
			int index = distribution(generator);
			indices_in_ascending_order.insert(index);
			if (indices_in_ascending_order.size() > size_set)
			{
				size_set = indices_in_ascending_order.size();
				counter += 1;
			}
		} while (counter < 3);


		// get pairs from queue
		std::vector<std::tuple<int, int, Scalar, Scalar>> triplet;
		int nb_of_elements_to_pop = *(indices_in_ascending_order.begin());
		auto it = indices_in_ascending_order.begin();
		counter = 0;
		while (counter < 3)
		{
			for (int i = 0; i < nb_of_elements_to_pop; i++)
				kpairs.pop();
			triplet.push_back(kpairs.top());
			nb_of_elements_to_pop = *std::next(it, 1) - *it;
			counter += 1;
			it = std::next(it, 1);
		}

		//return triplet
		RansacScheme::triplet t{ triplet[0], triplet[1], triplet[2] };
		vec_triplets.push_back(t);
	}
	return vec_triplets;
	
}


RansacScheme::triplet RansacScheme::pop_3_farthest_pairs()
{
	RansacScheme::triplet triplet;

	if (queue_ptr_->size() > 2)
	{
		triplet.pair1 = queue_ptr_->top();
		queue_ptr_->pop();
		
		find_2_farthest_pairs(triplet);
	}
	
	return triplet;
}


RansacScheme::triplet RansacScheme::pop_3_farthest_pairs(Scalar bbox_diag_ratio)
{
	RansacScheme::triplet triplet;

	if (queue_ptr_->size() > 2)
	{
		triplet.pair1 = queue_ptr_->top();
		queue_ptr_->pop();

		find_2_farthest_pairs_with_respect_to_geomVar(triplet, bbox_diag_ratio);
	}

	return triplet;
}

RansacScheme::triplet RansacScheme::pop_3_farthest_pairs(Scalar bbox_diag, Scalar lambda)
{
	RansacScheme::triplet triplet;

	if (queue_ptr_->size() > 2)
	{
		triplet.pair1 = queue_ptr_->top();
		queue_ptr_->pop();

		find_2_farthest_pairs_with_respect_to_geomVar(triplet, bbox_diag, lambda);
	}

	return triplet;
}


RansacScheme::pqueue RansacScheme::pop_k_farthest_pairs(RansacScheme::pqueue& queue, Scalar bbox_diag, int k, Scalar lambda)
{
	RansacScheme::pqueue queue_k_pairs;

	if ((int)queue.size() > k)
	{
		auto queue_copy = queue;
		std::tuple<int, int, Scalar, Scalar> pair1 = queue.top();
		queue.pop();
		
		// get k initial indices
		std::multimap<int, std::pair<Scalar, std::tuple<int, Scalar, Scalar>>> k_indices_and_sumdist;
		k_indices_and_sumdist.insert(std::make_pair(std::get<0>(pair1), std::make_pair(0.0, std::make_tuple(std::get<1>(pair1), std::get<2>(pair1), std::get<3>(pair1)))));
		while((int) k_indices_and_sumdist.size() < k)
		{
			auto pair = queue.top();
			queue.pop();
			k_indices_and_sumdist.insert(std::make_pair(std::get<0>(pair), std::make_pair(0.0, std::make_tuple(std::get<1>(pair), std::get<2>(pair), std::get<3>(pair)))));
		}

		// compute k-points metric (linear combination of cost and L2 distance)
		for (auto it_i = k_indices_and_sumdist.begin(); it_i != k_indices_and_sumdist.end(); it_i++)
		{
			Point first_k_point = std::get<0>(target_pointMap_->find_index(it_i->first)->second);
			Scalar priority_first_k_point = std::get<2>(target_pointMap_->find_index(it_i->first)->second);
			Scalar sum = 0.0;
			for (auto it_j = k_indices_and_sumdist.begin(); it_j != k_indices_and_sumdist.end(); it_j++)
			{
				Point second_k_point = std::get<0>(target_pointMap_->find_index(it_j->first)->second);
				Scalar priority_second_k_point = std::get<2>(target_pointMap_->find_index(it_j->first)->second);
				Scalar cost = priority_first_k_point * priority_second_k_point;
				sum += sqrt(compute_points_dist(first_k_point, second_k_point) / bbox_diag) + lambda*((max_cost_ - cost)/(max_cost_ - min_cost_));
			}
			
			it_i->second.first = sum;
		}

		// update k-points
			// find farther points
		while (!queue.empty())
		{
			auto pair = queue.top();
			queue.pop();
			// get point query
			Point& query_point = std::get<0>(target_pointMap_->find_index(std::get<0>(pair))->second);
			Scalar priority_query_point = std::get<2>(target_pointMap_->find_index(std::get<0>(pair))->second);
			// compute distance to the closest k-point from query point
			std::vector<std::pair<Scalar, int>> dists(k, std::make_pair(0.0, 0));
			int count = 0;
			for (auto it_ = k_indices_and_sumdist.begin(); it_ != k_indices_and_sumdist.end(); it_++)
			{
				Point& k_point = std::get<0>(target_pointMap_->find_index(it_->first)->second);
				int index_k_point = target_pointMap_->find_index(it_->first)->first;
				Scalar priority_k_point = std::get<2>(target_pointMap_->find_index(index_k_point)->second);
				Scalar cost = priority_k_point * priority_query_point;
				//std::cout << "cost : " << cost << std::endl;
				//std::cout << " dist : " << compute_points_dist(query_point, k_point) << std::endl;
				dists[count] = std::make_pair(sqrt(compute_points_dist(query_point, k_point)/bbox_diag) + lambda*((max_cost_ - cost)/(max_cost_ - min_cost_)), index_k_point);
				count += 1;
				//std::cout << "normalized cost : " << lambda*((max_cost_ - cost) / (max_cost_ - min_cost_)) << std::endl;
				//std::cout << "normalized dist : " << compute_points_dist(query_point, k_point) / bbox_diag << std::endl;
				//std::cout << " priority k point : " << priority_k_point << " , priority query_point : " << priority_query_point << " , cost : " << priority_k_point*priority_query_point << std::endl;
			}

			std::pair<Scalar, int> min_index = compute_min(dists, std::get<0>(pair1));
			// update k-points set (if needed)
			Point& point_to_be_replaced = std::get<0>(target_pointMap_->find_index(dists[min_index.second].second)->second);
			Scalar priority_point_to_be_replaced = std::get<2>(target_pointMap_->find_index(dists[min_index.second].second)->second);
			Scalar sum = 0.0;
			for (int i = 0; i < (int)dists.size(); i++)
				sum += dists[i].first;
			sum -= dists[min_index.second].first;

			if (sum > k_indices_and_sumdist.find(dists[min_index.second].second)->second.first)
			{
				for (auto it_ = k_indices_and_sumdist.begin(); it_ != k_indices_and_sumdist.end(); it_++)
				{
					Point& k_point = std::get<0>(target_pointMap_->find_index(it_->first)->second);
					Scalar priority_k_point = std::get<2>(target_pointMap_->find_index(it_->first)->second);
					Scalar cost_query = priority_k_point*priority_query_point;
					Scalar cost_to_replace = priority_point_to_be_replaced*priority_k_point;
					it_->second.first = it_->second.first + sqrt(compute_points_dist(query_point, k_point)/bbox_diag) + lambda*((max_cost_ - cost_query)/(max_cost_ - min_cost_)) - compute_points_dist(point_to_be_replaced, k_point)/bbox_diag - lambda*((max_cost_ - cost_to_replace)/(max_cost_ - min_cost_));
				}

				k_indices_and_sumdist.erase(k_indices_and_sumdist.find(dists[min_index.second].second));
				k_indices_and_sumdist.insert(std::make_pair(std::get<0>(pair), std::make_pair(sum, std::make_tuple(std::get<1>(pair), std::get<2>(pair), std::get<3>(pair)))));
			}
				
		}
	
		// fill vector of k farthest pair
		for (auto it_ = k_indices_and_sumdist.begin(); it_ != k_indices_and_sumdist.end(); it_++)
			queue_k_pairs.push(std::make_tuple(it_->first, std::get<0>(it_->second.second), std::get<1>(it_->second.second), std::get<2>(it_->second.second)));

		// update priority queue
		auto queue_k_pairs_copy = queue_k_pairs;
		while (!queue_copy.empty())
		{
			if (!queue_k_pairs_copy.empty())
			{
				auto pair = queue_copy.top();
				auto k_pair = queue_k_pairs_copy.top();
				queue_copy.pop();
				if (std::get<0>(pair) == std::get<0>(k_pair) && std::get<1>(pair) == std::get<1>(k_pair) && std::get<2>(pair) == std::get<2>(k_pair) && std::get<3>(pair) == std::get<3>(k_pair))
					queue_k_pairs_copy.pop();
				else
					queue.push(pair);
			}
			else
			{
				auto pair = queue_copy.top();
				queue_copy.pop();
				queue.push(pair);
			}
		}
	}

	return queue_k_pairs;
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

	// compute scale (target/source)
	Scalar scale = compute_scale(t);

	Eigen::Matrix3d scaling_mat = Eigen::Matrix3d::Zero();
	scaling_mat(0, 0) = 1.0 / scale;
	scaling_mat(1, 1) = 1.0 / scale;
	scaling_mat(2, 2) = 1.0 / scale;

	// points 
	VectorType target_pos1 = std::get<0>(target_pointMap_->find_index(std::get<0>(t.pair1))->second).pos();
	VectorType source_pos1 = std::get<0>(source_pointMap_->find_index(std::get<1>(t.pair1))->second).pos();
	VectorType target_pos2 = std::get<0>(target_pointMap_->find_index(std::get<0>(t.pair2))->second).pos();
	VectorType source_pos2 = std::get<0>(source_pointMap_->find_index(std::get<1>(t.pair2))->second).pos();
	VectorType target_pos3 = std::get<0>(target_pointMap_->find_index(std::get<0>(t.pair3))->second).pos();
	VectorType source_pos3 = std::get<0>(source_pointMap_->find_index(std::get<1>(t.pair3))->second).pos();

	// target and source data
	Eigen::Matrix3d source_data, target_data;
	source_data.col(0) = Eigen::Map<Eigen::ArrayXd>(source_pos1.data(), 3);
	source_data.col(1) = Eigen::Map<Eigen::ArrayXd>(source_pos2.data(), 3);
	source_data.col(2) = Eigen::Map<Eigen::ArrayXd>(source_pos3.data(), 3);
	target_data.col(0) = Eigen::Map<Eigen::ArrayXd>(target_pos1.data(), 3);
	target_data.col(1) = Eigen::Map<Eigen::ArrayXd>(target_pos2.data(), 3);
	target_data.col(2) = Eigen::Map<Eigen::ArrayXd>(target_pos3.data(), 3);

	// rescale target data
	target_data = scaling_mat*target_data;

	//computing centroids of both datasets
	Eigen::Vector3d centroid_source, centroid_target;
	centroid_source << (source_pos1(0,0) + source_pos2(0,0) + source_pos3(0,0)) / 3.0, (source_pos1(1,0) + source_pos2(1,0) + source_pos3(1,0)) / 3.0, (source_pos1(2,0) + source_pos2(2,0) + source_pos3(2,0)) / 3.0;
	centroid_target << (target_pos1(0,0) + target_pos2(0,0) + target_pos3(0,0)) / 3.0, (target_pos1(1,0) + target_pos2(1,0) + target_pos3(1,0)) / 3.0, (target_pos1(2,0) + target_pos2(2,0) + target_pos3(2,0)) / 3.0;

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
	translation = source_data.col(0) - (R*target_data).col(0); 

	// global transform 
	Eigen::Matrix4d transform = Eigen::Matrix4d::Zero(4, 4);
	transform.block(0, 0, 3, 3) = R*scaling_mat;
	transform.block(0, 3, 3, 1) = translation;
	transform(3, 3) = 1.0;

	//std::cout << "transform : " << transform << std::endl;

	return transform;
}


/* not sure this is the proper way to solve registration problem for affine transformation with 4 matching points. 
The method used (by solving a linear system) returns a 4-by-4 matrix (with coefficients for rotation, translation and scaling) but solves the 12 unknown as if they were decorrelated 
In fact, it is not quite the case. Some of the terms of the matrix are related to each other and we have only 10 unknowns [1 theta, 3 rot axis, 3 translation, 3 scaling ](and not 12)*/
Eigen::Matrix4d RansacScheme::compute_rigid_transform(RansacScheme::triplet t, std::tuple<int, int, Scalar, Scalar> q)
{
	// 4 pairs of points 
	VectorType target_pos1 = std::get<0>(target_pointMap_->find_index(std::get<0>(t.pair1))->second).pos();
	VectorType source_pos1 = std::get<0>(source_pointMap_->find_index(std::get<1>(t.pair1))->second).pos();
	VectorType target_pos2 = std::get<0>(target_pointMap_->find_index(std::get<0>(t.pair2))->second).pos();
	VectorType source_pos2 = std::get<0>(source_pointMap_->find_index(std::get<1>(t.pair2))->second).pos();
	VectorType target_pos3 = std::get<0>(target_pointMap_->find_index(std::get<0>(t.pair3))->second).pos();
	VectorType source_pos3 = std::get<0>(source_pointMap_->find_index(std::get<1>(t.pair3))->second).pos();
	VectorType target_q_pos = std::get<0>(target_pointMap_->find_index(std::get<0>(q))->second).pos();
	VectorType source_q_pos = std::get<0>(source_pointMap_->find_index(std::get<1>(q))->second).pos();

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
	Eigen::Matrix4d transform_transpose = (target_data.transpose()).colPivHouseholderQr().solve(source_data.transpose());


	//std::cout << "transform : " << transform_transpose.transpose() << std::endl;

	return transform_transpose.transpose();
}


Scalar RansacScheme::registrationErr(Eigen::Matrix4d& transform, RansacScheme::triplet& triplet)
{
	Eigen::Matrix3d R = transform.block(0, 0, 3, 3);
	Eigen::Vector3d T = transform.block(0, 3, 3, 1);
	Scalar err1 = ((Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair1))->second).pos().data(), 3, 1) - R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair1))->second).pos().data(), 3, 1) - T).cwiseAbs()).sum();
	Scalar err2 = ((Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair2))->second).pos().data(), 3, 1) - R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair2))->second).pos().data(), 3, 1) - T).cwiseAbs()).sum();
	Scalar err3 = ((Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair3))->second).pos().data(), 3, 1) - R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair3))->second).pos().data(), 3, 1) - T).cwiseAbs()).sum();

	//std::cout << "reg err : " << (err1 + err2 + err3) / 3.0;
	return (err1 + err2 + err3) / 3.0;
}

Scalar RansacScheme::registrationErr(Eigen::Matrix4d& transform, RansacScheme::triplet& triplet, std::tuple<int, int, Scalar, Scalar>& q)
{
	Eigen::Matrix3d R = transform.block(0, 0, 3, 3);
	Eigen::Vector3d T = transform.block(0, 3, 3, 1);
	Scalar err1 = ((Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair1))->second).pos().data(), 3, 1) - R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair1))->second).pos().data(), 3, 1) - T).cwiseAbs()).sum();
	Scalar err2 = ((Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair2))->second).pos().data(), 3, 1) - R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair2))->second).pos().data(), 3, 1) - T).cwiseAbs()).sum();
	Scalar err3 = ((Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair3))->second).pos().data(), 3, 1) - R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair3))->second).pos().data(), 3, 1) - T).cwiseAbs()).sum();
	Scalar errq = ((Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(q))->second).pos().data(), 3, 1) - R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(q))->second).pos().data(), 3, 1) - T).cwiseAbs()).sum();
	
	//std::cout << "reg err : " << (err1 + err2 + err3 + errq) / 4.0;
	return (err1 + err2 + err3 + errq) / 4.0;
}

Scalar RansacScheme::compute_angle(VectorType v1, VectorType v2)
{
	Scalar dot_product = v1.dot(v2);
	Scalar norm_v1 = v1.norm();
	Scalar norm_v2 = v2.norm();
	return acos(dot_product / (norm_v1*norm_v2)); // en radians
}

Scalar RansacScheme::normalErr(Eigen::Matrix4d& transform, RansacScheme::triplet& triplet)
{

	Eigen::Matrix3d R = transform.block(0, 0, 3, 3);
	Scalar err1 = compute_angle(Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair1))->second).normal().data(), 3, 1) ,  R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair1))->second).normal().data(), 3, 1));
	Scalar err2 = compute_angle(Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair2))->second).normal().data(), 3, 1), R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair2))->second).normal().data(), 3, 1));
	Scalar err3 = compute_angle(Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair3))->second).normal().data(), 3, 1), R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair3))->second).normal().data(), 3, 1));

	//std::cout << "  normal err : " << (err1 + err2 + err3) / 3.0 << std::endl;
	return (err1 + err2 + err3) / 3.0;
}

Scalar RansacScheme::normalErr(Eigen::Matrix4d& transform, RansacScheme::triplet& triplet, std::tuple<int, int, Scalar, Scalar>& q)
{

	Eigen::Matrix3d R = transform.block(0, 0, 3, 3);
	Scalar err1 = compute_angle(Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair1))->second).normal().data(), 3, 1), R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair1))->second).normal().data(), 3, 1));
	Scalar err2 = compute_angle(Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair2))->second).normal().data(), 3, 1), R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair2))->second).normal().data(), 3, 1));
	Scalar err3 = compute_angle(Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair3))->second).normal().data(), 3, 1), R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair3))->second).normal().data(), 3, 1));
	Scalar errq = compute_angle(Eigen::Map<Eigen::MatrixXd>(std::get<0>(source_pointMap_->find_index(std::get<1>(q))->second).normal().data(), 3, 1), R*Eigen::Map<Eigen::MatrixXd>(std::get<0>(target_pointMap_->find_index(std::get<0>(q))->second).normal().data(), 3, 1));

	//std::cout << "  normal err : " << (err1 + err2 + err3 + errq) / 4.0 << std::endl;
	return (err1 + err2 + err3 + errq) / 4.0;
}


bool RansacScheme::is_q_unique(RansacScheme::triplet t, std::tuple<int, int, Scalar, Scalar> q)
{
	/*std::cout << "target point 1 : " << std::get<0>(t.pair1).pos().transpose() << " , target point 2 : " << std::get<0>(t.pair2).pos().transpose() << " , target point 3 : " << std::get<0>(t.pair3).pos().transpose() << std::endl;
	std::cout << "source point 1 : " << std::get<1>(t.pair1).pos().transpose() << " , source point 2 : " << std::get<1>(t.pair2).pos().transpose() << " , source point 3 : " << std::get<1>(t.pair3).pos().transpose() << std::endl;
	std::cout << "target q : " << std::get<0>(q).pos().transpose() << " , source q : " << std::get<1>(q).pos().transpose() << std::endl;*/

	if (std::get<0>(target_pointMap_->find_index(std::get<0>(t.pair1))->second).pos() == std::get<0>(target_pointMap_->find_index(std::get<0>(q))->second).pos() || std::get<0>(target_pointMap_->find_index(std::get<0>(t.pair2))->second).pos() == std::get<0>(target_pointMap_->find_index(std::get<0>(q))->second).pos() || std::get<0>(target_pointMap_->find_index(std::get<0>(t.pair3))->second).pos() == std::get<0>(target_pointMap_->find_index(std::get<0>(q))->second).pos())
		return false;
	if (std::get<0>(source_pointMap_->find_index(std::get<1>(t.pair1))->second).pos() == std::get<0>(source_pointMap_->find_index(std::get<1>(q))->second).pos() || std::get<0>(source_pointMap_->find_index(std::get<1>(t.pair2))->second).pos() == std::get<0>(source_pointMap_->find_index(std::get<1>(q))->second).pos() || std::get<0>(source_pointMap_->find_index(std::get<1>(t.pair3))->second).pos() == std::get<0>(source_pointMap_->find_index(std::get<1>(q))->second).pos())
		return false;
	return true;
}


bool RansacScheme::is_valid(std::tuple<int, int, Scalar, Scalar> q, RansacScheme::triplet t, Scalar max_err_reg, Scalar max_err_norm)
{
	if (!is_q_unique(t, q))
		return false;
	Eigen::Matrix4d transform = compute_rigid_transform(t, q);
	Scalar reg_err = registrationErr(transform, t, q);
	Scalar norm_err = normalErr(transform, t, q);
	if (!(reg_err < max_err_reg && norm_err < max_err_norm))
		return false;
	else
		return true;

}


Scalar RansacScheme::compute_scale(RansacScheme::triplet& triplet)
{
	Point ps1, ps2, ps3, pt1, pt2, pt3;
	Scalar a, b, c, A, B, C;
	pt1 = std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair1))->second);
	ps1 = std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair1))->second);
	pt2 = std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair2))->second);
	ps2 = std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair2))->second);
	pt3 = std::get<0>(target_pointMap_->find_index(std::get<0>(triplet.pair3))->second);
	ps3 = std::get<0>(source_pointMap_->find_index(std::get<1>(triplet.pair3))->second);

	a = compute_points_dist(pt1, pt2);
	b = compute_points_dist(pt2, pt3);
	c = compute_points_dist(pt1, pt3);
	A = compute_points_dist(ps1, ps2);
	B = compute_points_dist(ps2, ps3);
	C = compute_points_dist(ps1, ps3);

	Scalar scale = (a/A + b/B + c/C) / 3.0;
	return scale;
}


Eigen::Matrix4d RansacScheme::ransac_algorithm(int nb_iterations, Scalar max_err_scale, Scalar max_err_reg, Scalar max_err_norm, Scalar bbox_diag_ratio, Scalar lambda, int k, std::string debug_filename)
{
	int counter = 0;
	while (counter < nb_iterations && queue_ptr_->size() > k /*2*/)
	{
		//triplet t = pop_triplet();
		//triplet t = pop_3_farthest_pairs(bbox_diag_ratio);
		int index_first_point = std::get<0>(queue_ptr_->top());
		pqueue kpairs_queue =  pop_k_farthest_pairs(*queue_ptr_, bbox_diag_ratio, k, lambda);
		if (counter ==0)
			write_ply_file_for_debug(debug_filename, target_pointMap_, kpairs_queue, index_first_point);
		std::vector<triplet> vec_triplets = pick_triplets_between_k_pairs(kpairs_queue, k, 50);
		for (int i = 0; i < 50; i++)
		{
			triplet t = vec_triplets[i];
			//triplet t = pick_triplet(8, 9, 23);
			/*std::cout << "triplet : " << std::endl;
			std::cout << std::get<0>(target_pointMap_->find_index(std::get<0>(t.pair1))->second).pos().transpose() << " , " << std::get<0>(source_pointMap_->find_index(std::get<1>(t.pair1))->second).pos().transpose() << " , " << std::get<2>(t.pair1) << std::endl;
			std::cout << std::get<0>(target_pointMap_->find_index(std::get<0>(t.pair2))->second).pos().transpose() << " , " << std::get<0>(source_pointMap_->find_index(std::get<1>(t.pair2))->second).pos().transpose() << " , " << std::get<2>(t.pair2) << std::endl;
			std::cout << std::get<0>(target_pointMap_->find_index(std::get<0>(t.pair3))->second).pos().transpose() << " , " << std::get<0>(source_pointMap_->find_index(std::get<1>(t.pair3))->second).pos().transpose() << " , " << std::get<2>(t.pair3) << std::endl;*/
			Scalar err_scale = scaleDiff(t).first;
			//std::cout << "err scale : " << err_scale << std::endl;
			if (err_scale < max_err_scale)
			{
				//Scalar avgScale = scaleDiff(t).second;
				Eigen::Matrix4d M = compute_rigid_transform(t);

				if (registrationErr(M, t) < max_err_reg && normalErr(M, t) < max_err_norm)
				{
					queue_copy_ = *queue_ptr_;
					while (!queue_copy_.empty())
					{
						std::tuple<int, int, Scalar, Scalar> q = queue_copy_.top();
						queue_copy_.pop();
						if (is_valid(q, t, max_err_reg, max_err_norm))
							return compute_rigid_transform(t, q);
					}
				}

			}
			counter += 1;
		}
	}
	return Eigen::Matrix4d::Identity();
}