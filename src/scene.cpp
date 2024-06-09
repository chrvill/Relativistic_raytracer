#include <Eigen/Dense>
#include <cmath>
#include <omp.h>
#include <cfenv>
#include "Image.h"
#include "metric.h"
#include "colorCalculator.h"
#include "scene.h"

Eigen::MatrixXd Scene::calculate_initial_rays(const Eigen::Vector3d& camera_dir, const Eigen::Vector3d& up_vector) {
    Eigen::Vector3d normalized_camera_dir = camera_dir.normalized();
    Eigen::Vector3d normalized_up_vector = up_vector.normalized();

    Eigen::Vector3d right_vector = normalized_camera_dir.cross(normalized_up_vector).normalized();

    double image_plane_height = 2 * focal_length * std::tan(fov / 2);
    double image_plane_width = aspect_ratio * image_plane_height;

    Eigen::MatrixXd ray_directions = Eigen::MatrixXd::Zero(image_width * image_height, 3);

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

img::ImageRGBf Scene::simulate_camera_rays(int n_steps, double h0, size_t image_height, size_t image_width) {
    rays = Eigen::MatrixXd::Zero(num_pixels, 8);

    double r     = camera_pos(0);
    double theta = camera_pos(1);
    double phi   = camera_pos(2);

    img::ImageRGBf image(image_height, image_width);
    //Eigen::MatrixXd colors = Eigen::MatrixXd::Zero(num_pixels, 3);

    double max_value = 0.0;

    #pragma omp parallel for
    for (int j = 0; j < num_pixels; ++j)
    {   
        Vector8d y0;
        Eigen::Vector3d ray = metric.transform_cartesian_vec(initial_rays.row(j), r, theta, phi);

        y0 << 0, r, theta, phi, 0, ray(0), ray(1), ray(2);
    
        double p0 = metric.compute_p0(y0(1), y0(2), y0(5), y0(6), y0(7));

        y0(4) = p0;
        double affine_parameter = 0.0;

        bool outside_celestial_sphere = false;
        bool below_EH = false;
        bool inside_disk = false;
        double redshift = 1.0;

        Vector8d y = y0;

        double h = h0;

        double optical_depth = 0.0;
        Eigen::Vector3d old_pos = metric.pos_to_cartesian(y(1), y(2), y(3));

        Eigen::Vector3d color(0.0, 0.0, 0.0);

        double T = 5000.0;
        for (int i = 0; i < n_steps; ++i) {
            Eigen::Vector3d new_pos = old_pos;

            y = metric.RKF45(y, h);

            affine_parameter += h;
            
            bool break_integration = metric.break_integration(y, outside_celestial_sphere, below_EH);

            if (render_disk) {
                if (y(1) <= 10.0 and h > 0.1) {
                    h = 0.1;
                }

                double z = y(1)*std::cos(y(2));
                inside_disk = std::abs(z) <= 10.0 and y(1) >= 5.0 and y(1) <= 10.0;

                if (inside_disk) {
                    Eigen::Vector3d new_pos = metric.pos_to_cartesian(y(1), y(2), y(3));
                    double dl = (new_pos - old_pos).norm();

                    double rho = compute_density(y(1), y(2), y(3));
                    optical_depth += dl*rho;

                    Eigen::Vector4d u_observer;
                    double u1_observer = 0.0;
                    double u2_observer = 0.0;
                    double u3_observer = metric.orbital_velocity(y(1));

                    double u0_observer = metric.compute_p0(y(1), y(2), u1_observer, u2_observer, u3_observer, -1);
                    u_observer << u0_observer, u1_observer, u2_observer, u3_observer;

                    // We're computing the redshift of a photon emitted by the camera and received at the disk
                    // so need to take the reciprocal to find the redshift of the photon emitted by the disk and received at the camera
                    redshift = 1.0/metric.compute_redshift(y0, y, u_observer, camera_u);

                    color += colorCalculator.compute_blackbody_RGB(T/redshift)*dl*std::exp(-optical_depth);
                }

                //double z = y(1)*std::cos(y(2));
                //inside_disk = (std::abs(z) <= 0.1) and (y(1) >= 5.0) and (y(1) <= 10.0);
            }
            if (break_integration or optical_depth >= 1.0) {
                break;
            }
        }

        /*
        Vector8d result = metric.solve_geodesic(y0, n_steps, h, affine_parameter, \
                                                outside_celestial_sphere, below_EH, inside_disk);
        */

        if (outside_celestial_sphere) 
        {
            Eigen::Vector3d cartesian_ray = metric.transform_vec_to_cartesian(Eigen::Vector3d(y(5), y(6), y(7)), y(1), y(2), y(3));

            double v_phi = std::atan2(cartesian_ray(1), cartesian_ray(0));
            double v_theta = std::acos(cartesian_ray(2)/cartesian_ray.norm());

            Eigen::Vector2d uv;;

            uv(0) = std::fmod(v_theta, M_PI)/M_PI;
            uv(1) = std::fmod(v_phi, 2*M_PI)/(2*M_PI);

            color = lookup_background(uv);
        }

        //if (color.maxCoeff() > 1) {
            //        color /= color.maxCoeff();
            //}

        for (int channel; channel < 3; ++channel) {
            if (color(channel) < 0)
            {
                color(channel) = 0;
            }
        }

        int index_x = j % image_width;
        int index_y = j / image_width;
        image(index_y, index_x) << color(0), color(1), color(2);

        //max_value = std::max({max_value, color(0), color(1), color(2)});
    }

    /*
    for (int y = 0; y < image_height; ++y) {
        for (int x = 0; x < image_width; ++x) {
            image(y, x) /= max_value;
        }
    }
    */

    img::save("renders/blackhole.png", image);
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
