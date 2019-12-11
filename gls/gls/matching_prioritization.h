#pragma once

#include <vector>
#include <string>
#include <tuple>
#include <queue>
#include "base.h"

typedef Point::Scalar Scalar;
class mycomparison
{
	bool reverse;
public:
	mycomparison(const bool& revparam = false);
	bool operator() (const std::tuple<Point, Point, Scalar, Scalar>& t1, const std::tuple<Point, Point, Scalar, Scalar>& t2) const;
};

typedef Point::Scalar Scalar;
class pair_priority_queue
{
private:
	std::priority_queue<std::tuple<Point, Point, Scalar, Scalar>, std::vector<std::tuple<Point, Point, Scalar, Scalar>>, mycomparison> queue_;

public:
	pair_priority_queue() {};
	pair_priority_queue(const std::tuple<Point, Point, Scalar, Scalar>& pair);
	pair_priority_queue(const std::vector<std::tuple<Point, Point, Scalar, Scalar>>& vector_pair);
	void add_pair(const std::tuple<Point, Point, Scalar, Scalar>& pair);
	void empty_queue();
};

typedef Point::Scalar Scalar;
std::pair<Point, Scalar> compute_point_priority(std::map<Point, std::vector<Scalar>, PointComp>::iterator point_geom_var, int nb_samples, Scalar alpha);

typedef Point::Scalar Scalar;
std::tuple<Point, Point, Scalar> compute_pair_cost(std::map<Point, std::vector<Scalar>, PointComp>::iterator source_geom_var, std::map<Point, std::vector<Scalar>, PointComp>::iterator target_geom_var, int nb_samples, Scalar alpha);