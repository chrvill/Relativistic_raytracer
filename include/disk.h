#pragma once 

#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include "metric.h"

class Disk {
public:
    Metric& metric;

    std::vector<Eigen::Vector3d> voronoi_points;
    int n_voronoi_points;

    double r_inner_edge;
    double r_outer_edge;

    Disk(Metric& metric_, int n_points, double rmin, double rmax) : metric(metric_) {
        r_inner_edge = rmin;
        r_outer_edge = rmax;

        n_voronoi_points = n_points;
        generate_voronoi_points(n_points);
    }

    void generate_voronoi_points(int n_points);

    Eigen::Vector3d get2DRotation(const Eigen::Vector3d& A, double angle, double angular_displacement);

    std::vector<double> getDiskDensity(const Eigen::Vector3d& sampleOrig, double r, double a, const Eigen::Vector3d& PosMin, double rmin);

    bool inline inside_disk(double r, double theta) {
        double z = r * std::cos(theta);
        return ((std::abs(z) <= 0.5) and (r >= r_inner_edge and r <= r_outer_edge));
    }
};
