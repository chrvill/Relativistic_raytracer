#pragma once 

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <omp.h>
#include <string>
#include <random>
#include "Image.h"
#include "metric.h"
#include "colorCalculator.h"
#include "disk.h"

class Scene {
public:
    double focal_length;
    size_t image_width;
    size_t image_height;
    double fov;
    int num_pixels;
    double aspect_ratio;

    Eigen::MatrixXd initial_rays;
    Eigen::Vector3d camera_pos;
    Metric& metric;
    Disk& disk;
    Eigen::MatrixXd rays;
    Eigen::Vector4d camera_u;
    ColorCalculator colorCalculator;

    bool render_disk;

    img::ImageRGBf background;
    std::string background_image_filename;

    Scene(double focal_length_, size_t image_width_, size_t image_height_, double fov_, Metric& metric_, Disk& disk_, std::string& cie_filename, std::string& background_image_filename_, bool render_disk_ = true): 
    metric(metric_), disk(disk_), colorCalculator(cie_filename) {
        focal_length = focal_length_;
        image_width = image_width_;
        image_height = image_height_;
        fov = fov_;
        num_pixels = image_width * image_height;
        aspect_ratio = static_cast<double>(image_width)/image_height;
        background_image_filename = background_image_filename_;
        render_disk = render_disk_;

        bool loaded = img::load(background_image_filename, background);
    }

    // Computing the directions of the rays that will be casted from the camera.
    Eigen::MatrixXd calculate_initial_rays(const Eigen::Vector3d& camera_dir, const Eigen::Vector3d& up_vector);

    void initialize(const Eigen::Vector3d camera_v, const Eigen::Vector3d& camera_pos_, 
                    const Eigen::Vector3d& camera_dir, const Eigen::Vector3d& up_vector);

    // Solve the geodesic equation for each of the rays, compute redshift and color for each pixel.
    img::ImageRGBf simulate_camera_rays(int n_steps, double h, size_t image_height, size_t image_width);

    Eigen::Vector3d lookup_background(Eigen::Vector2d uv);

    double compute_density(const double r, const double theta, const double phi);
};