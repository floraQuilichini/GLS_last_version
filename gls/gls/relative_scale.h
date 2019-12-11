#pragma once
# include <vector>
# include "base.h"

typedef Point::Scalar Scalar;
Scalar Median(std::vector<Scalar>::iterator begin, std::vector<Scalar>::iterator end);

typedef Point::Scalar Scalar;
typedef Point::VectorType VectorType;
Point find_point_matching(std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>> & p_profile, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& pts_profiles, VectorType w, Scalar alpha, Scalar nb_samples);

typedef Point::Scalar Scalar;
typedef Point::VectorType VectorType;
std::vector<std::pair<Point, Point>> compute_matching_pairs_after_rescaling(std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& source_gls_profiles, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& target_gls_profiles, VectorType w, Scalar alpha, Scalar nb_samples, bool cross_check=false);

typedef Point::Scalar Scalar;
typedef Point::VectorType VectorType;
std::pair<Scalar, std::vector<std::tuple<Point, Point, Scalar>>> compute_3_closest_pairs_with_scale(std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& source_gls_profiles, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& target_gls_profiles, Scalar ratio, int nb_samples, VectorType w, Scalar alpha = 1.0);

typedef Point::Scalar Scalar;
typedef Point::VectorType VectorType;
std::pair<Scalar, std::vector<std::tuple<Point, Point, Scalar>>> relative_scale_estimation(std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& source_gls_profiles, std::vector<std::pair<Point, std::vector<std::tuple<Scalar, Scalar, Scalar>>>>& target_gls_profiles, Scalar ratio, int nb_samples, VectorType w = VectorType::Ones(), Scalar alpha = 1.0);