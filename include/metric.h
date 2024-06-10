#pragma once

#include <Eigen/Dense>
#include <cmath>

using Vector8d = Eigen::Matrix<double, 8, 1>;

class Metric {
/*
General metric class that can be inherited by other metric classes.

The metric class contains the basic functions that are needed to compute the geodesic equations and solve them using the Runge-Kutta-Fehlberg method.

The default metric assumed here is the Minkowski metric in spherical coordinates, but this can be changed by inheriting from this class and overriding the functions '
g_tt, g_tph, g_rr, g_thth, g_phph and geodesic_eq_rhs. 

break_integration specifies the conditions under which the integration should be stopped. By default, it is set to false, but it can be overridden in the inherited class.

The coordinate transformation functions are just the standard transformations from spherical coordinates to cartesian and vice versa. These can also be overridden in the inherited class.
*/

public:
    double a = 0.0;
    double sign_a = 1.0;


    Metric(): a() {
        if (a != 0.0)
        {
            sign_a = a/std::abs(a);
        }
    } 

    // Components of the metric tensor
    virtual inline double g_tt(double r, double theta);

    virtual inline double g_tph(double r, double theta);

    virtual inline double g_rr(double r, double theta);

    virtual inline double g_thth(double r, double theta);

    virtual inline double g_phph(double r, double theta);

    // Computing the initial value of p^t = dx^t/dlambda. 
    // mu = p^\mu p_\mu
    double compute_p0(double r, double theta, double p1, double p2, double p3, double mu = 0);

    // Right-hand side of the geodesic equation
    virtual Vector8d geodesic_eq_rhs(const Vector8d& y);

    // Runge-Kutta-Fehlberg integration scheme
    virtual Vector8d RKF45(const Vector8d& y, double &h, double tol = 1e-5);

    // Conditions for breaking the integration
    virtual bool break_integration(const Vector8d& y, bool &outside_celestial_sphere, bool &below_EH);

    // Solve the geodesic equation for a given set of initial conditions
    Vector8d solve_geodesic(const Vector8d& y0, int n_steps, double h, double &affine_parameter, 
                            bool& outside_celestial_sphere, bool& below_EH, bool& inside_disk);

    // Take a position expressed in spherical/Boyer-Lindquist coordinates to Cartesian coordinates
    virtual Eigen::Vector3d pos_to_cartesian(double r, double theta, double phi);

    // Transformation matrix taking a vector expressed in Cartesian coordinates to spherical/Boyer-Lindquist coordinates
    virtual Eigen::Matrix3d transformationMatrix(double r, double theta, double phi);

    Eigen::Vector3d transform_cartesian_vec(const Eigen::Vector3d& vec, double r, double theta, double phi);

    Eigen::Vector3d transform_vec_to_cartesian(const Eigen::Vector3d& vec, double r, double theta, double phi);
    // Compute the redshift of a photon emitted at a position and with a momentum given by initial_y and received at a position and with a momentum given by final_y
    // The four-velocities of the emitter and receiver are given by u_emitter and u_observer, respectively
    double compute_redshift(const Vector8d& initial_y, const Vector8d& final_y, 
                            const Eigen::Vector4d& u_observer, const Eigen::Vector4d& u_emitter);

    // Analytical expression for the orbital velocity of a particle at a given radius (assuming equatorial orbit, and Schwarzschild/Kerr)
    double orbital_velocity(const double r);

    Eigen::Matrix4d lorentz_transformation(const Eigen::Vector3d& v);
};