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

bool mycomparison::operator()(const std::tuple<Point, Point, Scalar, Scalar>& t1, const std::tuple<Point, Point, Scalar, Scalar>& t2) const
{
	if (reverse) return (std::get<3>(t1)>std::get<3>(t2));
	else return (std::get<3>(t1)<std::get<3>(t2));
}


pair_priority_queue::pair_priority_queue(const std::tuple<Point, Point, Scalar, Scalar>& pair)
{
	queue_.push(pair);
}

pair_priority_queue::pair_priority_queue(const std::vector<std::tuple<Point, Point, Scalar, Scalar>>& vec_pair)
{
	for (int i=0; i<vec_pair.size(); i++)
		queue_.push(vec_pair[i]);
}

void pair_priority_queue::add_pair(const std::tuple<Point, Point, Scalar, Scalar>& pair)
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

std::pair<Point, Scalar> compute_point_priority(std::map<Point, std::vector<Scalar>, PointComp>::iterator point_geom_var, int nb_samples, Scalar alpha)
{
	Scalar sum = 0.0;
	for (int i = 0; i < nb_samples; i++)
		sum = sum + 1.0 - tanh(alpha*(point_geom_var->second)[i]);
	sum = sum / nb_samples;
	return std::make_pair(point_geom_var->first, sum);
}



std::tuple<Point, Point, Scalar> compute_pair_cost(std::map<Point, std::vector<Scalar>, PointComp>::iterator source_geom_var, std::map < Point, std::vector<Scalar>, PointComp > ::iterator target_geom_var, int nb_samples, Scalar alpha)
{
	std::pair<Point, Scalar> source_priority = compute_point_priority(source_geom_var, nb_samples, alpha);
	std::pair<Point, Scalar> target_priority = compute_point_priority(target_geom_var, nb_samples, alpha);

	return std::make_tuple(source_priority.first, target_priority.first, source_priority.second*target_priority.second);

}
