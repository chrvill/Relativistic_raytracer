#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <omp.h>
#include <cstdlib>
#include <string>

#include "kerr.h"
#include "colorCalculator.h"
#include "scene.h"
#include "Image.h"  // https://github.com/ThibaultLejemble/img
#include "disk.h"

int main() {
    // ------ Position of the camera ------
    double r     = 30.0; 
    double theta = M_PI / 2.0 - 0.1; 
    double phi   = M_PI/2; 
    double a     = 0.0; 

    // ------ Camera variables ------
    double focal_length = 1;
    constexpr size_t image_width = 1000;
    constexpr size_t image_height = 1000;
    constexpr size_t num_pixels = image_width * image_height;
    double fov = M_PI / 4;

    Eigen::Vector3d camera_pos(r, theta, phi);
    Eigen::Vector3d camera_dir(0.0, -1.0, 0.0);
    Eigen::Vector3d up_vector(0.0, 0.0, -1.0);
    Eigen::Vector3d camera_v(0.0, 0.0, 0.0);   // Local 3-velocity of the camera

    // ------ Kerr metric ------
    Kerr metric(0.999);

    // ------ Initializing the scene ------
    std::string cie_filename = "txtfiles/cie_interpolated.txt";
    std::string background_image_filename = "backgrounds/starmap_2020_8k.png";

    Disk disk(metric);

    bool render_disk = false;
    Scene scene(focal_length, image_width, image_height, fov, metric, disk, cie_filename, background_image_filename, render_disk);
    scene.initialize(camera_v, Eigen::Vector3d(r, theta, phi), camera_dir, up_vector);
    
    // ------ Simulating the camera rays and drawing the image ------
    img::ImageRGBf image = scene.simulate_camera_rays(1000, 1e-1, image_height, image_width);

    std::cout << "\nDone.\n";
}
