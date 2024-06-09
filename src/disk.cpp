#include <Eigen/Dense>
#include <iostream>
#include <cmath>

#include "disk.h"

double Disk::compute_density(const double r, const double theta, const double phi) {
    //double z = std::abs(r * std::cos(theta));

    //return 1e-1*(1.0/(r*std::sin(theta)) * std::exp(-5*z));

    Eigen::Vector3d cartesian = metric.pos_to_cartesian(r, theta, phi);

    return 1e3*perlin.noise(cartesian(0), cartesian(1), cartesian(2));
}

