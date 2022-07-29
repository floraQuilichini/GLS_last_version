#include "stdafx.h"

#include <cmath>
#include <algorithm>
#include <vector>
#include <tuple>
#include <string>
#include <queue> 
#include "Patate/grenaille.h"
#include "Eigen/Eigen"
#include "matching_prioritization.h"


mycomparison::mycomparison(const bool& revparam)
{
	reverse = revparam;
}

bool mycomparison::operator()(const std::tuple<int, int, Scalar, Scalar>& t1, const std::tuple<int, int, Scalar, Scalar>& t2) const
{
	if (reverse) return (std::get<3>(t1)<std::get<3>(t2));
	else return (std::get<3>(t1)>std::get<3>(t2));
}


pair_priority_queue::pair_priority_queue(const std::tuple<int, int, Scalar, Scalar>& pair)
{
	queue_.push(pair);
}

pair_priority_queue::pair_priority_queue(const std::vector<std::tuple<int, int, Scalar, Scalar>>& vec_pair)
{
	for (int i=0; i<(int)vec_pair.size(); i++)
		queue_.push(vec_pair[i]);
}

pair_priority_queue::pair_priority_queue(const std::vector<std::tuple<int, int, Scalar, Scalar>>& vec_pair, Scalar min_source_scale, Scalar min_target_scale, Scalar base)
{
	
	for (int i = 0; i < (int)vec_pair.size(); i++)
	{
		// compute scale 
		Scalar lag = std::get<2>(vec_pair[i]);
		Scalar scale = (min_target_scale / min_source_scale)*pow(base, -lag);
		// fill queue
		queue_.push(std::make_tuple(std::get<0>(vec_pair[i]), std::get<1>(vec_pair[i]), scale, std::get<3>(vec_pair[i])));
	}

}

pair_priority_queue::pair_priority_queue(const std::vector<std::tuple<int, int, Scalar, Scalar, Scalar>>& vec_pair, Scalar min_source_scale, Scalar min_target_scale, Scalar base)
{

	for (int i = 0; i < (int)vec_pair.size(); i++)
	{
		// compute scale 
		Scalar lag = std::get<2>(vec_pair[i]);
		Scalar scale = (min_target_scale / min_source_scale)*pow(base, -lag);

		// fill queue
		queue_.push(std::make_tuple(std::get<0>(vec_pair[i]), std::get<1>(vec_pair[i]), scale, std::get<3>(vec_pair[i])));
	}

}

pair_priority_queue::pair_priority_queue(const std::vector<std::tuple<std::pair<int, int>, Scalar, Scalar, Scalar>>& vec_pair, PointMap* source_pointMap, PointMap* target_pointMap, Scalar min_source_scale, Scalar min_target_scale, Scalar base)
{
	source_pointMap_ = source_pointMap;
	target_pointMap_ = target_pointMap;
	for (int i = 0; i < (int)vec_pair.size(); i++)
	{
		// compute scale 
		Scalar lag = std::get<1>(vec_pair[i]);
		Scalar scale = (min_target_scale / min_source_scale)*pow(base, -lag);

		// fill queue
		queue_.push(std::make_tuple(std::get<0>(vec_pair[i]).first, std::get<0>(vec_pair[i]).second, scale, std::get<2>(vec_pair[i])));
	}

}


void pair_priority_queue::add_pair(const std::tuple<int, int, Scalar, Scalar>& pair)
{
	queue_.push(pair);
}

void pair_priority_queue::empty_queue()
{
	while (!queue_.empty())
	{
		queue_.pop();
	}
}

void pair_priority_queue::display_queue()
{
	std::priority_queue<std::tuple<int, int, Scalar, Scalar>, std::vector<std::tuple<int, int, Scalar, Scalar>>, mycomparison> queue_copie = queue_;
	while (!queue_copie.empty())
	{
		std::tuple<int, int, Scalar, Scalar> tuple = queue_copie.top();
		queue_copie.pop();
		std::cout << "p_target : " << std::get<0>(target_pointMap_->find_index(std::get<0>(tuple))->second).pos() <<", p_source : " << std::get<0>(source_pointMap_->find_index(std::get<1>(tuple))->second).pos() << ", relative scale : " << std::get<2>(tuple) << ", cost : " << std::get<3>(tuple) << std::endl;
	}
}

size_t pair_priority_queue::get_queue_size() { return queue_.size(); }

std::priority_queue<std::tuple<int, int, Scalar, Scalar>, std::vector<std::tuple<int, int, Scalar, Scalar>>, mycomparison>* pair_priority_queue::get_queue_ptr()
{
	return (&queue_);
}

std::pair<Point, Scalar> compute_point_priority(std::map<Point, std::vector<Scalar>>::iterator point_geom_var, int nb_samples, Scalar alpha)
{
	Scalar sum = 0.0;
	for (int i = 0; i < nb_samples; i++)
		sum = sum + 1.0 - tanh(alpha*(point_geom_var->second)[i]);
	sum = sum / nb_samples;
	return std::make_pair(point_geom_var->first, sum);
}


Scalar compute_point_priority(std::vector<Scalar>& point_geom_var, int nb_samples, Scalar alpha)
{
	Scalar sum = 0.0;
	for (int i = 0; i < nb_samples; i++)
		sum = sum + 1.0 - tanh(alpha*(point_geom_var[i]));
	sum = sum / nb_samples;
	return sum;
}


std::tuple<Point, Point, Scalar> compute_pair_cost(std::map<Point, std::vector<Scalar>>::iterator source_geom_var, std::map < Point, std::vector<Scalar>> ::iterator target_geom_var, int nb_samples, Scalar alpha)
{
	std::pair<Point, Scalar> source_priority = compute_point_priority(source_geom_var, nb_samples, alpha);
	std::pair<Point, Scalar> target_priority = compute_point_priority(target_geom_var, nb_samples, alpha);

	return std::make_tuple(source_priority.first, target_priority.first, source_priority.second*target_priority.second);

}


Scalar compute_pair_cost(std::vector<Scalar>& source_geom_var, std::vector<Scalar>& target_geom_var, int nb_samples, Scalar alpha)
{
	Scalar source_priority = compute_point_priority(source_geom_var, nb_samples, alpha);
	Scalar target_priority = compute_point_priority(target_geom_var, nb_samples, alpha);

	return source_priority*target_priority;
}
