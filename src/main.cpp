#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <omp.h>
#include <memory>
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

std::unique_ptr<Metric> create_metric(const std::string& metric_name, const json& metric_parameters) {
    if (metric_name == "Kerr") {
        double a = metric_parameters["a"];
        return std::make_unique<Kerr>(a);
    }
    else if (metric_name == "Schwarzschild") {
        return std::make_unique<Schwarzschild>();
    }
    else if (metric_name == "Minkowski") {
        return std::make_unique<Metric>();
    }
    else {
        throw std::invalid_argument("Invalid metric name.");
    }
}

void create_image(json& scene_file) {
    auto camera = scene_file["camera"];
    auto cam_pos = camera["position"];

    auto metric_object = scene_file["metric"];
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
    std::string background_image_filename = scene_file["background"];

    auto metric_parameters = metric_object["parameters"];

    auto simulation_settings = scene_file["simulation_settings"];
    size_t n_steps = scene_file["simulation_settings"]["n_steps"];
    double h0 = scene_file["simulation_settings"]["initial_step_size"];
    bool render_disk = scene_file["render_disk"];
    double error_tolerance = simulation_settings["error_tolerance"];

    double r_inner_edge = scene_file["disk_parameters"]["radius_inner_edge"];
    double r_outer_edge = scene_file["disk_parameters"]["radius_outer_edge"];
    int n_voronoi_points = scene_file["disk_parameters"]["number_of_voronoi_points"];

    try {
        auto metric_pointer = create_metric(metric_name, metric_parameters);
        Metric& metric = *metric_pointer;

        Disk disk(metric, n_voronoi_points, r_inner_edge, r_outer_edge);

        Scene scene(focal_length, image_width, image_height, fov, metric, disk, cie_filename, background_image_filename, render_disk, error_tolerance);
        scene.initialize(camera_v, camera_pos, camera_dir, up_vector);

        // ------ Simulating the camera rays and drawing the image ------
        img::ImageRGBf image = scene.simulate_camera_rays(n_steps, h0);
        img::save(scene_file["output_file"], image);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
    }
}

int main() {
    omp_set_num_threads(16);

    auto start = std::chrono::high_resolution_clock::now();

    json scene_file = read_json("scenes/kerr.json");

    create_image(scene_file);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> render_time = end - start;
    std::cout << "Rendering time: " << render_time.count() << " s\n";

    std::cout << "\nDone.\n";

    return 0;
}
