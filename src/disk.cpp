#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include "PerlinNoise.h" //https://github.com/sol-prog/Perlin_Noise

#include "disk.h"

PerlinNoise noise;

// Sigmoid function
double sigmoid(double x) {
    return 1.0 / (1 + std::exp(-x));
}

double hscript_turb(const Eigen::Vector3d& p, int octaves, const PerlinNoise& noise) {
    double value = 0.0;
    double amplitude = 1.0;
    double frequency = 1.0;
    double max = 0.0;  // Used for normalizing result to [-1,1]

    for (int i = 0; i < octaves; i++) {
        value += noise.noise(p.x() * frequency, p.y() * frequency, p.z() * frequency) * amplitude;
        max += amplitude;
        amplitude *= 0.5;
        frequency *= 2.0;
    }

    return value / max;
}

// 2D Rotation function
Eigen::Vector3d Disk::get2DRotation(const Eigen::Vector3d& A, double angle, double angular_displacement) {
    double s = std::sin(angle + angular_displacement);
    double c = std::cos(angle + angular_displacement);
    Eigen::Matrix2d rotation_matrix;
    rotation_matrix << c, -s, s, c;
    Eigen::Vector2d A2 = rotation_matrix * Eigen::Vector2d(A.x(), A.y());
    return Eigen::Vector3d(A2.x(), A2.y(), A.z());
}

void Disk::generate_voronoi_points(int n_points) {
    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_r(r_inner_edge, r_outer_edge);
    std::uniform_real_distribution<> dis_phi(0, 2 * M_PI);

    // Generate points
    for (int i = 0; i < n_points; ++i) {
        double r = dis_r(gen);
        double phi = dis_phi(gen);
        
        // Convert spherical/Boyer-Lindquist to Cartesian coordinates
        Eigen::Vector3d pos = metric.pos_to_cartesian(r, M_PI / 2, phi);
        
        // Store the point in the voronoi_points vector
        voronoi_points.emplace_back(pos(0), pos(1), 0); // z = 0 as the points are in the xy-plane
    }
}

std::vector<double> Disk::getDiskDensity(const Eigen::Vector3d& point, double r, double a, const Eigen::Vector3d& closest_point, double shortest_distance) {
    // Field Falloffs
    double innerFallOff = 1.0 / (1 + std::exp(-1 * (r - r_inner_edge)));
    double outerFallOff = 1.0 / (1 + std::exp(0.3 * (r - r_outer_edge)));
    double detailFallOff = 1.0 / (1 + std::exp(0.5 * (r - r_outer_edge)));
    
    // Initialize Coordinates
    Eigen::Vector3d tC = get2DRotation(point, std::sqrt(1.0 / std::pow(r + 3, 2)) * -((a < 0) ? -1 : 1) * 96, 10);
    Eigen::Vector3d pC = (tC - closest_point) / 110.0;
    
    // Noise Field
    double f1 = 1.0 / (1 + shortest_distance);
    double element0 = std::abs(hscript_turb(pC*56, 5, noise));
    double element1 = std::abs(hscript_turb(pC.cwiseProduct(Eigen::Vector3d(12, 12, 56)), 5, noise));
    double element2 = std::pow(f1, 1.5) * f1;
    double element3 = f1 * 0.005 * detailFallOff;
    double element4 = std::pow(std::abs(hscript_turb((pC + Eigen::Vector3d(2, 55, 32)).cwiseProduct(Eigen::Vector3d(12, 12, 56)), 5, noise)), 2) * 7;
    double density = 5;
    double EmissionField = ((((element0 + element1) / 2.0) * element2 * density) + element3) * outerFallOff;
    double AbsorptionField = ((((element0 + element1 + element4) / 3.0) * element2 * density) + element3) * outerFallOff;

    // Break-up upper and lower edge
    Eigen::Vector3d adjustedSampleOrig = point - Eigen::Vector3d(0, 0, innerFallOff * ((point.z() < 0) ? -1 : 1));
    double sampleAngle = std::asin((adjustedSampleOrig.normalized()).dot(Eigen::Vector3d(0, 0, ((point.z() < 0) ? -1 : 1))));

    double diskSlope = M_PI / 50; 
    double slopeNoise = EmissionField * 5;
    diskSlope += slopeNoise;
    
    // In-Exact Bounds test
    if (sampleAngle < diskSlope) {
        // Calculate smooth vertical gradient
        double x = std::abs(point.z());
        double h = std::tan(diskSlope) * r;
        double A = 2; // Gradient Amplitude
        double N_0 = 2.7; // First Power
        double N_1 = 7; // Second Power
        double k_0 = std::pow(A, -N_0);
        double k_1 = std::pow(h, -N_1) * h;
        double m = -(A / h);

        if (std::abs(k_1 * std::pow(x, N_1)) < h) {
            double gx = std::pow(((std::abs(k_1 * std::pow(x, N_1)) - h) * m), N_0) * k_0 * A;
            
            // Absorption Falloff
            double absFalloff = 1 + (10.0 / (1 + std::exp(-0.2 * (r - r_outer_edge))));
            //std::cout << EmissionField * gx << "\n";
            return {EmissionField * gx, std::pow(AbsorptionField, 1.3) * gx * absFalloff * 3};
        } else {
            return {0.0, 0.0};
        }
    } else {
        return {0.0, 0.0};
    }
}