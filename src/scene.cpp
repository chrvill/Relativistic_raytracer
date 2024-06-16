#include <Eigen/Dense>
#include <cmath>
#include <omp.h>
#include <cfenv>
#include <mutex>
#include "Image.h"
#include "metric.h"
#include "disk.h"
#include "colorCalculator.h"
#include "scene.h"

void print_progress_bar(float progress) {
    int barWidth = 70;
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

Eigen::MatrixXd Scene::calculate_initial_rays(const Eigen::Vector3d& camera_dir, const Eigen::Vector3d& up_vector) {
    Eigen::Vector3d normalized_camera_dir = camera_dir.normalized();
    Eigen::Vector3d normalized_up_vector = up_vector.normalized();

    Eigen::Vector3d right_vector = normalized_camera_dir.cross(normalized_up_vector).normalized();

    double image_plane_height = 2 * focal_length * std::tan(fov / 2);
    double image_plane_width = aspect_ratio * image_plane_height;

    Eigen::MatrixXd ray_directions = Eigen::MatrixXd::Zero(image_width * image_height, 3);

    #pragma omp parallel for
    for (int y = 0; y < image_height; ++y) {
        for (int x = 0; x < image_width; ++x) {
            int index = y * image_width + x;

            double pixel_x = ((x + 0.5) / image_width - 0.5) * image_plane_width;
            double pixel_y = ((y + 0.5) / image_height - 0.5) * image_plane_height;

            Eigen::Vector3d image_plane_point(pixel_x, pixel_y, focal_length);

            Eigen::Vector3d ray_direction_unnormalized = pixel_x * right_vector +
                                            pixel_y * normalized_up_vector +
                                            focal_length * normalized_camera_dir;

            Eigen::Vector3d ray_direction = ray_direction_unnormalized.normalized();

            ray_directions.row(index) = ray_direction;
        }
    }

    return ray_directions;
}

void Scene::initialize(const Eigen::Vector3d camera_v, const Eigen::Vector3d& camera_pos_, 
                        const Eigen::Vector3d& camera_dir, const Eigen::Vector3d& up_vector) {
    camera_pos = camera_pos_;
    initial_rays = calculate_initial_rays(camera_dir, up_vector);

    double u0_emitter = metric.compute_p0(camera_pos(0), camera_pos(1), camera_v(0), camera_v(1), camera_v(2), -1);
    camera_u << u0_emitter, camera_v(0), camera_v(1), camera_v(2);
}

img::ImageRGBf Scene::simulate_camera_rays(int n_steps, double h0) {
    rays = Eigen::MatrixXd::Zero(num_pixels, 8);

    double r     = camera_pos(0);
    double theta = camera_pos(1);
    double phi   = camera_pos(2);

    img::ImageRGBf image(image_height, image_width);
    
    Eigen::Vector3d local_camera_velocity = metric.compute_local_cartesian_velocity(Eigen::Vector3d(camera_u(1), camera_u(2), camera_u(3)), r, theta);
    
    Eigen::Matrix4d lorentz = metric.lorentz_transformation(local_camera_velocity);

    int progress = 0;
    int update_frequency = num_pixels/100;

    std::mutex progress_mutex;

    double g_tt = metric.g_tt(r, theta);
    double g_rr = metric.g_rr(r, theta);
    double g_thth = metric.g_thth(r, theta);
    double g_phph = metric.g_phph(r, theta);
    double g_tph = metric.g_tph(r, theta);

    Eigen::Matrix3d transformation = metric.transformationMatrix(r, theta, phi);

    double max_color_value = 0.0;

    #pragma omp parallel
    {
        int local_progress = 0;
        #pragma omp for
        for (int j = 0; j < num_pixels; ++j)
        {   
            Vector8d y0;
            Eigen::Vector3d cartesian_ray = initial_rays.row(j);

            double ray0 = cartesian_ray.norm();;

            Eigen::Vector4d ray_prime = lorentz * Eigen::Vector4d(ray0, cartesian_ray(0), cartesian_ray(1), cartesian_ray(2));

            //Eigen::Vector3d ray_local = metric.transform_cartesian_vec(Eigen::Vector3d(ray_prime(1), ray_prime(2), ray_prime(3)), r, theta, phi);

            Eigen::Vector3d ray_local = transformation * Eigen::Vector3d(ray_prime(1), ray_prime(2), ray_prime(3));

            //Eigen::Vector4d ray = metric.transform_vec_to_global(Eigen::Vector4d(ray_prime(0), ray_local(0), ray_local(1), ray_local(2)), r, theta, phi);
            Eigen::Vector4d ray(ray_prime(0), ray_local(0)/std::sqrt(g_rr), ray_local(1)/std::sqrt(g_thth), ray_local(2)/std::sqrt(g_phph));

            y0 << 0, r, theta, phi, 0, ray(1), ray(2), ray(3);

            //double p0 = metric.compute_p0(y0(1), y0(2), y0(5), y0(6), y0(7));
            double p0 = metric.compute_p0(g_tt, g_rr, g_thth, g_phph, g_tph, y0(5), y0(6), y0(7));

            y0(4) = p0;

            double affine_parameter = 0.0;

            bool outside_celestial_sphere = false;
            bool below_EH = false;
            bool inside_disk = false;

            double redshift = 1.0;

            Vector8d y = y0;

            double h = h0;

            double optical_depth = 0.0;

            Eigen::Vector3d color(0.0, 0.0, 0.0);
            Eigen::Vector4d u_observer(0.0, 0.0, 0.0, 0.0);

            double T = 5000.0;
            for (int i = 0; i < n_steps; ++i) {
                Eigen::Vector3d old_pos = metric.pos_to_cartesian(y(1), y(2), y(3));

                y = metric.RKF45(y, h, error_tolerance);

                Eigen::Vector3d new_pos = metric.pos_to_cartesian(y(1), y(2), y(3));

                affine_parameter += h;
                
                bool break_integration = metric.break_integration(y, outside_celestial_sphere, below_EH);

                if (break_integration) {
                    break;
                }

                if (render_disk) {
                    if (y(1) <= disk.r_outer_edge and h > 0.1) {
                        h = 0.1;
                    }

                    inside_disk = disk.inside_disk(y(1), y(2));

                    if (inside_disk) {
                        double dl = (new_pos - old_pos).norm();

                        Eigen::Vector3d closest_point(0.0, 0.0, 0.0);
                        double min_distance = 1e10;

                        Eigen::Vector3d point = disk.get2DRotation(new_pos, std::sqrt(1.0 / std::pow(y(1) + 3, 2)) * -((metric.a < 0) ? -1 : 1) * 96, 10);
                        for (int k = 0; k < disk.n_voronoi_points; ++k) {
                            Eigen::Vector3d voronoi_point = disk.voronoi_points[k];
                            //double distance = (new_pos - voronoi_point).norm();
                            double distance = (point - voronoi_point).norm();

                            if (distance < min_distance) {
                                min_distance = distance;
                                closest_point = voronoi_point;
                            }
                        }

                        std::vector<double> densities = disk.getDiskDensity(new_pos, y(1), metric.a, closest_point, min_distance);
                        double emission = densities[0];
                        double absorption = densities[1];

                        double u1_observer = 0.0;
                        double u2_observer = 0.0;
                        double u3_observer = metric.orbital_velocity(y(1));

                        double u0_observer = metric.compute_p0(y(1), y(2), u1_observer, u2_observer, u3_observer, -1);
                        u_observer << u0_observer, u1_observer, u2_observer, u3_observer;

                        // We're computing the redshift of a photon emitted by the camera and received at the disk
                        // so need to take the reciprocal to find the redshift of the photon emitted by the disk and received at the camera
                        redshift = 1.0/metric.compute_redshift(y0, y, u_observer, camera_u);
                        
                        //color = colorCalculator.compute_blackbody_RGB(T/redshift)*(1 + densities[0])*std::pow(redshift, 5);
                        optical_depth += dl*absorption;
                        Eigen::Vector3d raw_color = colorCalculator.compute_blackbody_RGB(T/redshift)*std::pow(redshift, 5);
                        color += raw_color*emission*dl*std::exp(-optical_depth);

                        old_pos = new_pos;

                        if (optical_depth > 1) break;
                    }
                }
            }

            if (color.maxCoeff() > max_color_value) {
                max_color_value = color.maxCoeff();
            }

            if (outside_celestial_sphere) 
            {
                Eigen::Vector3d cartesian_ray = metric.transform_vec_to_cartesian(Eigen::Vector3d(y(5), y(6), y(7)), y(1), y(2), y(3));

                double v_phi = std::atan2(cartesian_ray(1), cartesian_ray(0));
                double v_theta = std::acos(cartesian_ray(2)/cartesian_ray.norm());

                Eigen::Vector2d uv;

                uv(0) = std::fmod(v_theta, M_PI)/M_PI;
                uv(1) = std::fmod(v_phi, 2*M_PI)/(2*M_PI);

                color += lookup_background(uv);
            }

            color = color.cwiseMax(0.0);

            int index_x = j % image_width;
            int index_y = j / image_width;
            image(index_y, index_x) << color(0), color(1), color(2);

            local_progress++;

            if (local_progress % update_frequency == 0) {
                std::lock_guard<std::mutex> lock(progress_mutex);
                progress += update_frequency;
                float progress_ratio = static_cast<float>(progress) / num_pixels;
                print_progress_bar(progress_ratio);
            }
        }
        std::lock_guard<std::mutex> lock(progress_mutex);
        progress += local_progress % update_frequency;
        float progress_ratio = static_cast<float>(progress) / num_pixels;
        print_progress_bar(progress_ratio);
    }

    // Ensure the progress bar shows 100% at the end
    print_progress_bar(1.0);
    std::cout << std::endl;

    #pragma omp parallel for
    for (int i = 0; i < num_pixels; ++i) {
        int index_x = i % image_width;
        int index_y = i / image_width;
        image(index_y, index_x) = image(index_y, index_x)/max_color_value;
    }

    return image;
}

Eigen::Vector3d Scene::lookup_background(Eigen::Vector2d uv)
{
    uv(0) = static_cast<int>(std::round(uv(0)*background.height()));
    uv(1) = static_cast<int>(std::round(uv(1)*background.width()));

    Eigen::Vector3d color; 

    for (int channel = 0; channel < 3; ++channel)
    {
        color(channel) = background(uv(0), uv(1))[channel];        
    }

    return color;
}

double Scene::compute_density(const double r, const double theta, const double phi) {
    double z = std::abs(r * std::cos(theta));

    return 1e-1*(1.0/(r*std::sin(theta)) * std::exp(-5*z));
}
