#pragma once 

#include <cmath>
#include <array>
#include <Eigen/Dense>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cnpy.h>
#include <unsupported/Eigen/Splines>

using Array1000d = Eigen::Array<double, 1000, 1>;

class ColorCalculator {
public:
    double c   = 2.99792458e8;     // Speed of light
    double k_B = 1.3806504e-23;    // Boltzmann's constant
    double h   = 6.62607015e-34;   // Planck's constant
    std::array<double, 1000> lambdas;
    std::array<double, 1000> x_bar;
    std::array<double, 1000> y_bar;
    std::array<double, 1000> z_bar;

    Array1000d lambdas_eigen;
    Array1000d x_bar_eigen;
    Array1000d y_bar_eigen;
    Array1000d z_bar_eigen;

    size_t lambdas_size = 1000;

    Eigen::Matrix3d XYZ_to_RGB;

    ColorCalculator(const std::string& cie_filename);

    // Read the CIE data, which are samples of the color matching functions, from a file
    void readFile(const std::string& filename);

    // Compute the blackbody distribution for a given temperature and the range of wavelengths
    inline Array1000d blackbody_distribution(double T, const Array1000d& lambdas) {
        return 1.0 / (lambdas.cube() * lambdas.square())
            * 1.0 / (exp((h * c) / (lambdas * k_B * T)) - 1.0);
    }   

    // Compute the XYZ values for a given temperature
    Eigen::Vector3d compute_blackbody_XYZ(double T);

    // Compute the RGB color for a given temperature
    Eigen::Vector3d compute_blackbody_RGB(double T);
};