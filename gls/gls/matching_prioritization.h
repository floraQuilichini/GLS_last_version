#pragma once

#include <vector>
#include <string>
#include <tuple>
#include <queue>
#include "base.h"
#include "relative_scale.h"


class mycomparison
{
	bool reverse;
public:
	mycomparison(const bool& revparam = false);
	bool operator() (const std::tuple<int, int, Scalar, Scalar>& t1, const std::tuple<int, int, Scalar, Scalar>& t2) const;
};


class pair_priority_queue
{
private:
	std::priority_queue<std::tuple<int, int, Scalar, Scalar>, std::vector<std::tuple<int, int, Scalar, Scalar>>, mycomparison> queue_;
	PointMap* target_pointMap_ = nullptr;
	PointMap* source_pointMap_ = nullptr;

public:
	pair_priority_queue() {};
	pair_priority_queue(const std::tuple<int, int, Scalar, Scalar>& pair);
	pair_priority_queue(const std::vector<std::tuple<int, int, Scalar, Scalar>>& vector_pair);
	pair_priority_queue(const std::vector<std::tuple<int, int, Scalar, Scalar>>& vec_pair, Scalar min_source_scale, Scalar min_target_scale, Scalar base);
	pair_priority_queue(const std::vector<std::tuple<int, int, Scalar, Scalar, Scalar>>& vec_pair, Scalar min_source_scale, Scalar min_target_scale, Scalar base);
	pair_priority_queue(const std::vector<std::tuple<std::pair<int, int>, Scalar, Scalar, Scalar>>& vec_pair, PointMap* source_pointMap, PointMap* target_pointMap, Scalar min_source_scale, Scalar min_target_scale, Scalar base);
	void add_pair(const std::tuple<int, int, Scalar, Scalar>& pair);
	void empty_queue();
	void display_queue();
	size_t get_queue_size();
	std::priority_queue<std::tuple<int, int, Scalar, Scalar>, std::vector<std::tuple<int, int, Scalar, Scalar>>, mycomparison>* get_queue_ptr();
};


std::pair<Point, Scalar> compute_point_priority(std::map<Point, std::vector<Scalar>>::iterator point_geom_var, int nb_samples, Scalar alpha);


Scalar compute_point_priority(std::vector<Scalar>& point_geom_var, int nb_samples, Scalar alpha);


std::tuple<Point, Point, Scalar> compute_pair_cost(std::map<Point, std::vector<Scalar>>::iterator source_geom_var, std::map<Point, std::vector<Scalar>>::iterator target_geom_var, int nb_samples, Scalar alpha);


Scalar compute_pair_cost(std::vector<Scalar>& source_geom_var, std::vector<Scalar>& target_geom_var, int nb_samples, Scalar alpha);