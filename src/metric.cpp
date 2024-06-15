#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include "metric.h"

double Metric::compute_p0(double r, double theta, double p1, double p2, double p3, double mu) {
    double g_tt   = this->g_tt(r, theta);
    double g_tph  = this->g_tph(r, theta);
    double g_rr   = this->g_rr(r, theta);
    double g_thth = this->g_thth(r, theta);
    double g_phph = this->g_phph(r, theta);

    double g_frac = g_tph / g_tt;

    if (g_tt < 0) {
        return -g_frac * p3 + std::sqrt((g_frac) * (g_frac) * p3 * p3 - 1 / g_tt * (g_rr * p1 * p1 + g_thth * p2 * p2 + g_phph * p3 * p3 - mu));
    } else {
        return -g_frac * p3 - std::sqrt((g_frac) * (g_frac) * p3 * p3 - 1 / g_tt * (g_rr * p1 * p1 + g_thth * p2 * p2 + g_phph * p3 * p3 - mu));
    }
}

double Metric::compute_p0(double gtt, double grr, double gthth, double gphph, double gtph, double p1, double p2, double p3, double mu) {
    double g_frac = gtph / gtt;

    if (gtt < 0) {
        return -g_frac * p3 + std::sqrt((g_frac) * (g_frac) * p3 * p3 - 1 / gtt * (grr * p1 * p1 + gthth * p2 * p2 + gphph * p3 * p3 - mu));
    } else {
        return -g_frac * p3 - std::sqrt((g_frac) * (g_frac) * p3 * p3 - 1 / gtt * (grr * p1 * p1 + gthth * p2 * p2 + gphph * p3 * p3 - mu));
    }
}


void Metric::geodesic_eq_rhs(const Vector8d& y, Vector8d& derivatives) {
    //Vector8d derivatives; 

    double r = y(1);
    double theta = y(2);
    double p0 = y(4);
    double p1 = y(5);
    double p2 = y(6);
    double p3 = y(7);

    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    derivatives << p0, 
                   p1, 
                   p2, 
                   p3, 
                   0.0, 
                   r*p2*p2 + r*sin_theta*sin_theta*p3*p3, 
                   -2.0/r*p1*p2 + sin_theta*cos_theta*p3*p3, 
                   -2.0/r*p1*p3 - 2.0*cos_theta/sin_theta*p2*p3;

    //return derivatives;
}

Vector8d Metric::RKF45(const Vector8d& y, double &h, double tol) {
    // Initialize the error to a large value
    double abs_error = 1e100;

    Vector8d k1;
    Vector8d k2;
    Vector8d k3;
    Vector8d k4;
    Vector8d k5;
    Vector8d k6;

    Vector8d error;

    while (abs_error > tol)
    {
        geodesic_eq_rhs(y, k1);
        k1 *= h;

        geodesic_eq_rhs(y + B(0, 0) * k1, k2);
        k2 *= h;

        geodesic_eq_rhs(y + B(1, 0) * k1 + B(1, 1) * k2, k3);
        k3 *= h;

        geodesic_eq_rhs(y + B(2, 0) * k1 + B(2, 1) * k2 + B(2, 2) * k3, k4);
        k4 *= h;

        geodesic_eq_rhs(y + B(3, 0) * k1 + B(3, 1) * k2 + B(3, 2) * k3 + B(3, 3) * k4, k5);
        k5 *= h;

        geodesic_eq_rhs(y + B(4, 0) * k1 + B(4, 1) * k2 + B(4, 2) * k3 + B(4, 3) * k4 + B(4, 4) * k5, k6);
        k6 *= h;


        //Vector8d error = -1.0/360.0 * k1 + 128.0/4275.0 * k3 + 2197.0/75240.0 * k4 - 1.0/50.0 * k5 - 2.0/55.0 * k6;
        error = CT(0) * k1 + CT(1) * k2 + CT(2) * k3 + CT(3) * k4 + CT(4) * k5 + CT(5) * k6;

        //abs_error = std::sqrt(error(0) * error(0) + error(1) * error(1) + error(2) * error(2) + error(3) * error(3) + error(4) * error(4) + error(5) * error(5) + error(6) * error(6) + error(7) * error(7));
        abs_error = error.norm();

        // Adjust the step size based on the error estimate and the tolerance
        h = 0.9*h*std::pow(tol/abs_error, 1.0/5.0);
    }

    //return y + 16.0 / 135.0 * k1 + 6656.0 / 12825.0 * k3 + 28651.0 / 56430.0 * k4 - 9.0/50.0 * k5 + 2.0/55.0 * k6;
    return y + CH(0) * k1 + CH(1) * k2 + CH(2) * k3 + CH(3) * k4 + CH(4) * k5 + CH(5) * k6;
}

Eigen::Matrix3d Metric::transformationMatrix(double r, double theta, double phi) {
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    Eigen::Matrix3d M;

    M(0, 0) = sin_theta*cos_phi;
    M(0, 1) = sin_theta*sin_phi;
    M(0, 2) = cos_theta;

    M(1, 0) = cos_theta*cos_phi;
    M(1, 1) = cos_theta*sin_phi;
    M(1, 2) = -sin_theta;

    M(2, 0) = -sin_phi;
    M(2, 1) = cos_phi;
    M(2, 2) = 0;

    return M;
}

Eigen::Vector3d Metric::transform_cartesian_vec(const Eigen::Vector3d& vec, double r, double theta, double phi) {
    return transformationMatrix(r, theta, phi) * vec;
}

Eigen::Vector3d Metric::transform_vec_to_cartesian(const Eigen::Vector3d& vec, double r, double theta, double phi) {
    Eigen::Matrix3d M = transformationMatrix(r, theta, phi);

    Eigen::Vector3d new_vec(vec(0), vec(1), vec(2));

    new_vec(0) *= std::sqrt(g_rr(r, theta));
    new_vec(1) *= std::sqrt(g_thth(r, theta));
    new_vec(2) *= std::sqrt(g_phph(r, theta));

    return M.transpose() * new_vec;
}

double Metric::compute_redshift(const Vector8d& initial_y, const Vector8d& final_y, 
                            const Eigen::Vector4d& u_observer, const Eigen::Vector4d& u_emitter)
{
    double r_init = initial_y(1);
    double theta_init = initial_y(2);
    double p0_init = initial_y(4);
    double p1_init = initial_y(5);
    double p2_init = initial_y(6);
    double p3_init = initial_y(7);

    double r_final = final_y(1);
    double theta_final = final_y(2);
    double p0_final = final_y(4);
    double p1_final = final_y(5);
    double p2_final = final_y(6);
    double p3_final = final_y(7);

    double lambda_init = -1.0/(this->g_tt(r_init, theta_init)*p0_init*u_emitter(0) \
                                + this->g_rr(r_init, theta_init)*p1_init*u_emitter(1) \
                                + this->g_thth(r_init, theta_init)*p2_init*u_emitter(2) \
                                + this->g_phph(r_init, theta_init)*p3_init*u_emitter(3)\
                                + this->g_tph(r_init, theta_init)*(p0_init*u_emitter(3) + p3_init*u_emitter(0)));

    double lambda_final = -1.0/(this->g_tt(r_final, theta_final)*p0_final*u_observer(0) \
                                + this->g_rr(r_final, theta_final)*p1_final*u_observer(1) \
                                + this->g_thth(r_final, theta_final)*p2_final*u_observer(2) \
                                + this->g_phph(r_final, theta_final)*p3_final*u_observer(3)\
                                + this->g_tph(r_final, theta_final)*(p0_final*u_observer(3) + p3_final*u_observer(0)));

    return lambda_final/lambda_init;
}

double Metric::orbital_velocity(const double r) {
    return sign_a*std::sqrt(r)/(r*std::sqrt(r*r - 3*r + 2*std::abs(a)*std::sqrt(r)));
}

Eigen::Matrix4d Metric::lorentz_transformation(const Eigen::Vector3d& v) {
    double gamma = 1.0/std::sqrt(1 - v.dot(v));
    double beta = v.norm();

    Eigen::Matrix4d L;

    if (beta > 0.0) {
        L(0, 0) = gamma;
        L(0, 1) = -gamma*v(0);
        L(0, 2) = -gamma*v(1);
        L(0, 3) = -gamma*v(2);

        L(1, 0) = -gamma*v(0);
        L(1, 1) = 1 + (gamma - 1)*(v(0)*v(0))/(beta*beta);
        L(1, 2) = (gamma - 1)*(v(0)*v(1))/(beta*beta);
        L(1, 3) = (gamma - 1)*(v(0)*v(2))/(beta*beta);

        L(2, 0) = -gamma*v(1);
        L(2, 1) = (gamma - 1)*(v(1)*v(0))/(beta*beta);
        L(2, 2) = 1 + (gamma - 1)*(v(1)*v(1))/(beta*beta);
        L(2, 3) = (gamma - 1)*(v(1)*v(2))/(beta*beta);

        L(3, 0) = -gamma*v(2);
        L(3, 1) = (gamma - 1)*(v(2)*v(0))/(beta*beta);
        L(3, 2) = (gamma - 1)*(v(2)*v(1))/(beta*beta);
        L(3, 3) = 1 + (gamma - 1)*(v(2)*v(2))/(beta*beta);
    }
    else {
        L = Eigen::Matrix4d::Identity();
    }
    return L;
}

Eigen::Vector3d Metric::compute_local_cartesian_velocity(const Eigen::Vector3d& v, double r, double theta) {
    double gamma = compute_p0(r, theta, v(0), v(1), v(2), -1);
    
    Eigen::Vector3d v_cartesian = transform_vec_to_cartesian(v, r, theta, 0.0);

    return v_cartesian/gamma;
}

Eigen::Vector4d Metric::transform_vec_to_global(const Eigen::Vector4d& vec, double r, double theta, double phi) {
    Eigen::Vector4d vec_global;

    vec_global(0) = vec(0);
    vec_global(1) = vec(1)/std::sqrt(g_rr(r, theta));
    vec_global(2) = vec(2)/std::sqrt(g_thth(r, theta));
    vec_global(3) = vec(3)/std::sqrt(g_phph(r, theta));

    return vec_global;
}
