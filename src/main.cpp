#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <omp.h>
#include <cstdlib>
#include <string>
#include <chrono>
#include "nlohmann/json.hpp" // https://github.com/nlohmann/json

#include "metric.h"
#include "schwarzschild.h"
#include "kerr.h"
#include "colorCalculator.h"
#include "scene.h"
#include "Image.h"  // https://github.com/ThibaultLejemble/img
#include "disk.h"

using json = nlohmann::json;

json read_json(const std::string& scene_filename) {
    FILE* scene_file = std::fopen(scene_filename.c_str(), "r");
    char buffer[1024];
    std::string file; 

    // Read the file line by line
    while (std::fgets(buffer, sizeof(buffer), scene_file)) {
        file += buffer;
    }

    json data = json::parse(file);

    return data;
}

int main() {
    json data = read_json("scenes/minkowski.json");

    auto camera = data["camera"];
    auto cam_pos = camera["position"];

    auto metric_object = data["metric"];
    std::string metric_name = metric_object["name"];

    // ------ Camera variables ------
    auto camera_properties = camera["properties"];

    double focal_length = camera_properties["focal_length"];
    double fov = camera_properties["fov"];

    size_t image_width = camera_properties["image_width"];
    size_t image_height = camera_properties["image_height"];

    Eigen::Vector3d camera_pos(cam_pos["r"], cam_pos["theta"], cam_pos["phi"]);
    Eigen::Vector3d camera_dir(camera["camera_dir"][0], camera["camera_dir"][1], camera["camera_dir"][2]);
    Eigen::Vector3d up_vector(camera["up_vector"][0], camera["up_vector"][1], camera["up_vector"][2]);
    Eigen::Vector3d camera_v(camera["camera_velocity"][0], camera["camera_velocity"][1], camera["camera_velocity"][2]);   // Local 3-velocity of the camera
    
    // ------ Initializing the scene ------
    std::string cie_filename = "txtfiles/cie_interpolated.txt";
    std::string background_image_filename = data["background"];

    auto metric_parameters = metric_object["parameters"];

    auto simulation_settings = data["simulation_settings"];
    size_t n_steps = simulation_settings["n_steps"];
    double h0 = simulation_settings["initial_step_size"];
    bool render_disk = data["render_disk"];

    auto start = std::chrono::high_resolution_clock::now();

    if (metric_name == "Kerr") {

        double a = metric_parameters["a"];
        Kerr metric(a);
        std::cout << a << "\n";
        Disk disk(metric);

        Scene scene(focal_length, image_width, image_height, fov, metric, disk, cie_filename, background_image_filename, render_disk);
        scene.initialize(camera_v, camera_pos, camera_dir, up_vector);

        // ------ Simulating the camera rays and drawing the image ------
        img::ImageRGBf image = scene.simulate_camera_rays(n_steps, h0, image_height, image_width);
        img::save(data["output_file"], image);
    }
    else if (metric_name == "Schwarzschild") {
        Schwarzschild metric;
        Disk disk(metric);

        Scene scene(focal_length, image_width, image_height, fov, metric, disk, cie_filename, background_image_filename, render_disk);
        scene.initialize(camera_v, camera_pos, camera_dir, up_vector);

        // ------ Simulating the camera rays and drawing the image ------
        img::ImageRGBf image = scene.simulate_camera_rays(n_steps, h0, image_height, image_width);
        img::save(data["output_file"], image);
    }
    else if (metric_name == "Minkowski") {
        Metric metric;
        Disk disk(metric);
        
        Scene scene(focal_length, image_width, image_height, fov, metric, disk, cie_filename, background_image_filename, render_disk);
        scene.initialize(camera_v, camera_pos, camera_dir, up_vector);

        // ------ Simulating the camera rays and drawing the image ------
        img::ImageRGBf image = scene.simulate_camera_rays(n_steps, h0, image_height, image_width);
        img::save(data["output_file"], image);
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> render_time = end - start;
    std::cout << "Rendering time: " << render_time.count() << " s\n";

    std::cout << "\nDone.\n";

    return 0;
}
