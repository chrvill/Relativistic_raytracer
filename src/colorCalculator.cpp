#include <cmath>
#include <Eigen/Dense>
#include <cstdio>
#include <iostream>
#include <vector>
#include "colorCalculator.h"

ColorCalculator::ColorCalculator(const std::string& cie_filename) {
    readFile(cie_filename);

    XYZ_to_RGB << 3.24156456, -1.53766524 , -0.49870224,
                    -0.96920119,  1.87588535,  0.04155324,
                    0.05562416, -0.20395525 ,  1.05685902;
}

void ColorCalculator::readFile(const std::string& filename) {
    // Open the file in binary read mode
    FILE* file = std::fopen(filename.c_str(), "r");

    // Check if the file was successfully opened
    if (!file) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
    }

    // Create a buffer to hold a line of the file
    char buffer[1024];
    std::vector<std::vector<double>> data;

    // Read the file line by line
    while (std::fgets(buffer, sizeof(buffer), file)) {
        std::istringstream lineStream(buffer);
        std::vector<double> row;
        double value;

        // Read each value in the line and only take the first 4 columns
        for (int i = 0; i < 4 && (lineStream >> value); ++i) {
            row.push_back(value);
        }

        if (!row.empty()) {
            data.push_back(row);
        }
    }

    // Close the file
    std::fclose(file);

    for (size_t i = 0; i < data.size(); ++i) {
        lambdas[i] = data[i][0];
        x_bar[i] = data[i][1];
        y_bar[i] = data[i][2];
        z_bar[i] = data[i][3];
    }

    lambdas_eigen = Eigen::Map<const Eigen::Array<double, 1000, 1>>(lambdas.data());
    x_bar_eigen = Eigen::Map<const Eigen::Array<double, 1000, 1>>(x_bar.data());
    y_bar_eigen = Eigen::Map<const Eigen::Array<double, 1000, 1>>(y_bar.data());
    z_bar_eigen = Eigen::Map<const Eigen::Array<double, 1000, 1>>(z_bar.data());
}

Eigen::Vector3d ColorCalculator::compute_blackbody_XYZ(double T) {
    Array1000d I = blackbody_distribution(T, lambdas_eigen*1e-9);

    double dlambda = lambdas[1] - lambdas[0];

    double X = (I * x_bar_eigen).sum() * dlambda;
    double Y = (I * y_bar_eigen).sum() * dlambda;
    double Z = (I * z_bar_eigen).sum() * dlambda;

    return Eigen::Vector3d(X, Y, Z);
}

Eigen::Vector3d ColorCalculator::compute_blackbody_RGB(double T) {
    Eigen::Vector3d XYZ = compute_blackbody_XYZ(T);
    return XYZ_to_RGB*XYZ;
}
