#pragma once 

#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include "perlin.h"
#include "metric.h"

class Disk {
public:
    PerlinNoise perlin;
    Metric& metric;

    Disk(Metric& metric_) : metric(metric_) {
        perlin = PerlinNoise();
    }

    double compute_density(const double r, const double theta, const double phi);
};
