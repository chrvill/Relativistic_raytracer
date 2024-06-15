#include <Eigen/Dense>
#include <cmath>
#include "kerr.h"

inline double Kerr::Sigma(double r, double theta) {
    return r * r + a * a * std::cos(theta) * std::cos(theta);
}

inline double Kerr::Delta(double r, double theta) {
    return r * r - 2 * r + a * a;
}

inline double Kerr::Lambda(double r, double theta) {
    return (r * r + a * a) * (r * r + a * a) - a * a * Delta(r, theta) * std::sin(theta) * std::sin(theta);

}

inline double Kerr::g_tt(double r, double theta) {
    double sigma = r * r + a * a * std::cos(theta) * std::cos(theta);
    return -(1 - 2 * r / sigma);
}

inline double Kerr::g_tph(double r, double theta) {
    double sigma = r * r + a * a * std::cos(theta) * std::cos(theta);
    return -2 * a * r / sigma * std::sin(theta) * std::sin(theta);
}

inline double Kerr::g_rr(double r, double theta) {
    double sigma = r * r + a * a * std::cos(theta) * std::cos(theta);
    double delta = r * r - 2 * r + a * a;
    return sigma / delta;
}

inline double Kerr::g_thth(double r, double theta) {
    double sigma = r * r + a * a * std::cos(theta) * std::cos(theta);
    return sigma;
}

inline double Kerr::g_phph(double r, double theta) {
    double sigma = r * r + a * a * std::cos(theta) * std::cos(theta);
    double delta = r * r - 2 * r + a * a;
    double Lambda = (r * r + a * a) * (r * r + a * a) - a * a * delta * std::sin(theta) * std::sin(theta);
    return Lambda / sigma * std::sin(theta) * std::sin(theta);
}

void Kerr::geodesic_eq_rhs(const Vector8d& y, Vector8d& derivatives) {
    //Eigen::Matrix<double, 8, 1> derivatives;
    double t = y(0);
    double r = y(1);
    double theta = y(2);
    double phi = y(3);
    double p0 = y(4);
    double p1 = y(5);
    double p2 = y(6);
    double p3 = y(7);

    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);

    double cos_theta2 = cos_theta*cos_theta;
    double sin_theta2 = sin_theta*sin_theta;
    double sin_theta3 = sin_theta2*sin_theta;
    double sin_theta4 = sin_theta2*sin_theta2;

    double sin_cos = sin_theta*cos_theta;
    double sin2_cos2 = sin_theta2*cos_theta2;

    double a2 = a*a; 
    double a3 = a2*a;
    double a4 = a2*a2;

    double r2 = r*r;
    double r3 = r2*r;
    double r4 = r2*r2;

    double Sigma = r2 + a2*cos_theta2;
    double Sigma2 = Sigma*Sigma;
    double Sigma3 = Sigma2*Sigma;

    double Delta = -2*r + a2 + r2;
    double Delta2 = Delta*Delta;

    double ar = a*r;
    double ar2 = a*r2;

    double a2r2 = a2*r2;

    derivatives << p0, 
                   p1,
                   p2,
                   p3,
                   -2*p0*p1*(-1.0*ar*(4*ar2*sin_theta2/Sigma2 - 2*a*sin_theta2/(Sigma))/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4) + 0.5*(-4*r2/Sigma2 + 2/(Sigma))*(-2*a2*r*sin_theta2 - a4*cos_theta2 - a2r2*cos_theta2 - a2r2 - r4)/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4)) - 2*p0*p2*(2.0*a2*r*(-2*a2*r*sin_theta2 - a4*cos_theta2 - a2r2*cos_theta2 - a2r2 - r4)*sin_cos/(Sigma2*(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4)) - 1.0*ar*(-4*a3*r*sin_theta3*cos_theta/Sigma2 - 4*ar*sin_cos/(Sigma))/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4)) - 2*p1*p3*(-1.0*ar*(-4*a2r2*sin_theta2/Sigma2 + 2*a2*sin_theta2/(Sigma) + 2*r)*sin_theta2/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4) + 0.5*(4*ar2*sin_theta2/Sigma2 - 2*a*sin_theta2/(Sigma))*(-2*a2*r*sin_theta2 - a4*cos_theta2 - a2r2*cos_theta2 - a2r2 - r4)/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4)) - 2*p2*p3*(-1.0*ar*((4*a4*r*sin_theta3*cos_theta/Sigma2 + 4*a2*r*sin_cos/(Sigma))*sin_theta2 + 2*(2*a2*r*sin_theta2/(Sigma) + a2 + r2)*sin_cos)/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4) + 0.5*(-4*a3*r*sin_theta3*cos_theta/Sigma2 - 4*ar*sin_cos/(Sigma))*(-2*a2*r*sin_theta2 - a4*cos_theta2 - a2r2*cos_theta2 - a2r2 - r4)/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4)),
                   2.0*a2*p1*p2*sin_cos/(Sigma) + 1.0*r*p2*p2*Delta/(Sigma) - 0.5*p0*p0*(4*r2/Sigma2 - 2/(Sigma))*Delta/(Sigma) - 1.0*p0*p3*(-4*ar2*sin_theta2/Sigma2 + 2*a*sin_theta2/(Sigma))*Delta/(Sigma) - 0.5*p1*p1*(2*r/Delta + (2 - 2*r)*(Sigma)/Delta2)*Delta/(Sigma) + 0.5*p3*p3*Delta*(-4*a2r2*sin_theta2/Sigma2 + 2*a2*sin_theta2/(Sigma) + 2*r)*sin_theta2/(Sigma),
                   2.0*a2*r*p0*p0*sin_cos/Sigma3 - 1.0*a2*p1*p1*sin_cos/((Sigma)*Delta) + 1.0*a2*p2*p2*sin_cos/(Sigma) - 2.0*r*p1*p2/(Sigma) - 1.0*p0*p3*(4*a3*r*sin_theta3*cos_theta/Sigma2 + 4*ar*sin_cos/(Sigma))/(Sigma) - 0.5*p3*p3*(-(4*a4*r*sin_theta3*cos_theta/Sigma2 + 4*a2*r*sin_cos/(Sigma))*sin_theta2 - 2*(2*a2*r*sin_theta2/(Sigma) + a2 + r2)*sin_cos)/(Sigma),
                   -2*p0*p1*(-1.0*ar*(-4*r2/Sigma2 + 2/(Sigma))/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4) + 0.5*(4*ar2*sin_theta2/Sigma2 - 2*a*sin_theta2/(Sigma))*(-2*r + Sigma)/(2*a2*r*sin_theta4 - 2*a2*r*sin_theta2 - 2*r3*sin_theta2 + a4*sin2_cos2 + a2r2*sin2_cos2 + a2r2*sin_theta2 + r4*sin_theta2)) - 2*p0*p2*(-4.0*a3*r2*sin_cos/(Sigma2*(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4)) + 0.5*(-4*a3*r*sin_theta3*cos_theta/Sigma2 - 4*ar*sin_cos/(Sigma))*(-2*r + Sigma)/(2*a2*r*sin_theta4 - 2*a2*r*sin_theta2 - 2*r3*sin_theta2 + a4*sin2_cos2 + a2r2*sin2_cos2 + a2r2*sin_theta2 + r4*sin_theta2)) - 2*p1*p3*(-1.0*ar*(4*ar2*sin_theta2/Sigma2 - 2*a*sin_theta2/(Sigma))/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4) + 0.5*(-2*r + Sigma)*(-4*a2r2*sin_theta2/Sigma2 + 2*a2*sin_theta2/(Sigma) + 2*r)*sin_theta2/(2*a2*r*sin_theta4 - 2*a2*r*sin_theta2 - 2*r3*sin_theta2 + a4*sin2_cos2 + a2r2*sin2_cos2 + a2r2*sin_theta2 + r4*sin_theta2)) - 2*p2*p3*(-1.0*ar*(-4*a3*r*sin_theta3*cos_theta/Sigma2 - 4*ar*sin_cos/(Sigma))/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4) + 0.5*((4*a4*r*sin_theta3*cos_theta/Sigma2 + 4*a2*r*sin_cos/(Sigma))*sin_theta2 + 2*(2*a2*r*sin_theta2/(Sigma) + a2 + r2)*sin_cos)*(-2*r + Sigma)/(2*a2*r*sin_theta4 - 2*a2*r*sin_theta2 - 2*r3*sin_theta2 + a4*sin2_cos2 + a2r2*sin2_cos2 + a2r2*sin_theta2 + r4*sin_theta2));

    //derivatives(0) = p0;
    //derivatives(1) = p1;
    //derivatives(2) = p2;
    //derivatives(3) = p3;

    //derivatives(4) = -2*p0*p1*(-1.0*ar*(4*ar2*sin_theta2/Sigma2 - 2*a*sin_theta2/(Sigma))/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4) + 0.5*(-4*r2/Sigma2 + 2/(Sigma))*(-2*a2*r*sin_theta2 - a4*cos_theta2 - a2r2*cos_theta2 - a2r2 - r4)/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4)) - 2*p0*p2*(2.0*a2*r*(-2*a2*r*sin_theta2 - a4*cos_theta2 - a2r2*cos_theta2 - a2r2 - r4)*sin_cos/(Sigma2*(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4)) - 1.0*ar*(-4*a3*r*sin_theta3*cos_theta/Sigma2 - 4*ar*sin_cos/(Sigma))/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4)) - 2*p1*p3*(-1.0*ar*(-4*a2r2*sin_theta2/Sigma2 + 2*a2*sin_theta2/(Sigma) + 2*r)*sin_theta2/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4) + 0.5*(4*ar2*sin_theta2/Sigma2 - 2*a*sin_theta2/(Sigma))*(-2*a2*r*sin_theta2 - a4*cos_theta2 - a2r2*cos_theta2 - a2r2 - r4)/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4)) - 2*p2*p3*(-1.0*ar*((4*a4*r*sin_theta3*cos_theta/Sigma2 + 4*a2*r*sin_cos/(Sigma))*sin_theta2 + 2*(2*a2*r*sin_theta2/(Sigma) + a2 + r2)*sin_cos)/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4) + 0.5*(-4*a3*r*sin_theta3*cos_theta/Sigma2 - 4*ar*sin_cos/(Sigma))*(-2*a2*r*sin_theta2 - a4*cos_theta2 - a2r2*cos_theta2 - a2r2 - r4)/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4));

    //derivatives(5) = 2.0*a2*p1*p2*sin_cos/(Sigma) + 1.0*r*p2*p2*Delta/(Sigma) - 0.5*p0*p0*(4*r2/Sigma2 - 2/(Sigma))*Delta/(Sigma) - 1.0*p0*p3*(-4*ar2*sin_theta2/Sigma2 + 2*a*sin_theta2/(Sigma))*Delta/(Sigma) - 0.5*p1*p1*(2*r/Delta + (2 - 2*r)*(Sigma)/Delta2)*Delta/(Sigma) + 0.5*p3*p3*Delta*(-4*a2r2*sin_theta2/Sigma2 + 2*a2*sin_theta2/(Sigma) + 2*r)*sin_theta2/(Sigma);

    //derivatives(6) = 2.0*a2*r*p0*p0*sin_cos/Sigma3 - 1.0*a2*p1*p1*sin_cos/((Sigma)*Delta) + 1.0*a2*p2*p2*sin_cos/(Sigma) - 2.0*r*p1*p2/(Sigma) - 1.0*p0*p3*(4*a3*r*sin_theta3*cos_theta/Sigma2 + 4*ar*sin_cos/(Sigma))/(Sigma) - 0.5*p3*p3*(-(4*a4*r*sin_theta3*cos_theta/Sigma2 + 4*a2*r*sin_cos/(Sigma))*sin_theta2 - 2*(2*a2*r*sin_theta2/(Sigma) + a2 + r2)*sin_cos)/(Sigma);

    //derivatives(7) = -2*p0*p1*(-1.0*ar*(-4*r2/Sigma2 + 2/(Sigma))/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4) + 0.5*(4*ar2*sin_theta2/Sigma2 - 2*a*sin_theta2/(Sigma))*(-2*r + Sigma)/(2*a2*r*sin_theta4 - 2*a2*r*sin_theta2 - 2*r3*sin_theta2 + a4*sin2_cos2 + a2r2*sin2_cos2 + a2r2*sin_theta2 + r4*sin_theta2)) - 2*p0*p2*(-4.0*a3*r2*sin_cos/(Sigma2*(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4)) + 0.5*(-4*a3*r*sin_theta3*cos_theta/Sigma2 - 4*ar*sin_cos/(Sigma))*(-2*r + Sigma)/(2*a2*r*sin_theta4 - 2*a2*r*sin_theta2 - 2*r3*sin_theta2 + a4*sin2_cos2 + a2r2*sin2_cos2 + a2r2*sin_theta2 + r4*sin_theta2)) - 2*p1*p3*(-1.0*ar*(4*ar2*sin_theta2/Sigma2 - 2*a*sin_theta2/(Sigma))/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4) + 0.5*(-2*r + Sigma)*(-4*a2r2*sin_theta2/Sigma2 + 2*a2*sin_theta2/(Sigma) + 2*r)*sin_theta2/(2*a2*r*sin_theta4 - 2*a2*r*sin_theta2 - 2*r3*sin_theta2 + a4*sin2_cos2 + a2r2*sin2_cos2 + a2r2*sin_theta2 + r4*sin_theta2)) - 2*p2*p3*(-1.0*ar*(-4*a3*r*sin_theta3*cos_theta/Sigma2 - 4*ar*sin_cos/(Sigma))/(2*a2*r*sin_theta2 - 2*a2*r - 2*r3 + a4*cos_theta2 + a2r2*cos_theta2 + a2r2 + r4) + 0.5*((4*a4*r*sin_theta3*cos_theta/Sigma2 + 4*a2*r*sin_cos/(Sigma))*sin_theta2 + 2*(2*a2*r*sin_theta2/(Sigma) + a2 + r2)*sin_cos)*(-2*r + Sigma)/(2*a2*r*sin_theta4 - 2*a2*r*sin_theta2 - 2*r3*sin_theta2 + a4*sin2_cos2 + a2r2*sin2_cos2 + a2r2*sin_theta2 + r4*sin_theta2));

    //return derivatives;
}

Eigen::Vector3d Kerr::CartesianToBLVector(const Eigen::Vector3d& vector, double r, double theta, double phi, double a) {
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    double sqrt_factor = 1 / std::sqrt(r * r + a * a * cos_theta * cos_theta);

    Eigen::Matrix3d transform;
    transform(0, 0) = r * sin_theta * cos_phi * sqrt_factor;
    transform(0, 1) = r * sin_theta * sin_phi * sqrt_factor;
    transform(0, 2) = std::sqrt(r * r + a * a) * cos_theta * sqrt_factor;
    
    transform(1, 0) = std::sqrt(r * r + a * a) * cos_theta * cos_phi * sqrt_factor;
    transform(1, 1) = std::sqrt(r * r + a * a) * cos_theta * sin_phi * sqrt_factor;
    transform(1, 2) = -r * sin_theta * sqrt_factor;
    
    transform(2, 0) = -sin_phi;
    transform(2, 1) = cos_phi;
    transform(2, 2) = 0;

    Eigen::Vector3d BL_vector = transform * vector;

    BL_vector(0) *= 1.0/std::sqrt(g_rr(r, theta));
    BL_vector(1) *= 1.0/std::sqrt(g_thth(r, theta));
    BL_vector(2) *= 1.0/std::sqrt(g_phph(r, theta));

    return BL_vector;
}

Eigen::Vector3d Kerr::pos_to_cartesian(double r, double theta, double phi) {
    double x = std::sqrt(r*r + a*a) * std::sin(theta) * std::cos(phi);
    double y = std::sqrt(r*r + a*a) * std::sin(theta) * std::sin(phi);
    double z = r * std::cos(theta);

    return Eigen::Vector3d(x, y, z);
}

Eigen::Matrix3d Kerr::transformationMatrix(double r, double theta, double phi) {
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    double sqrt_factor = 1 / std::sqrt(r * r + a * a * cos_theta * cos_theta);

    Eigen::Matrix3d M;

    M(0, 0) = r*sqrt_factor*sin_theta*cos_phi;
    M(0, 1) = r*sqrt_factor*sin_theta*sin_phi;
    M(0, 2) = std::sqrt(r*r + a*a)*sqrt_factor*cos_theta;

    M(1, 0) = std::sqrt(r*r + a*a)*sqrt_factor*cos_theta*cos_phi;
    M(1, 1) = std::sqrt(r*r + a*a)*sqrt_factor*cos_theta*sin_phi;
    M(1, 2) = -r*sqrt_factor*sin_theta;

    M(2, 0) = -sin_phi;
    M(2, 1) = cos_phi;
    M(2, 2) = 0;

    return M;
}

Eigen::Vector3d Kerr::compute_local_cartesian_velocity(const Eigen::Vector3d& v, double r, double theta) {    
    double v0 = compute_p0(r, theta, v(0), v(1), v(2), -1);
    std::cout << "v0\t" << v0 << std::endl;

    double sigma = Sigma(r, theta);
    double delta = Delta(r, theta);
    double lambda = Lambda(r, theta);

    double u0 = std::sqrt(delta*sigma/lambda)*v0;
    double u1 = std::sqrt(sigma/delta)*v(0);
    double u2 = std::sqrt(sigma)*v(1);
    double u3 = -2.0*a*r*std::sin(theta)/std::sqrt(lambda*sigma)*v0 + std::sin(theta)*std::sqrt(lambda/sigma)*v(2);

    Eigen::Vector3d v_BL(u1/u0, u2/u0, u3/u0);

    return transform_vec_to_cartesian(v_BL, r, theta, 0.0);
}

Eigen::Vector4d Kerr::transform_vec_to_global(const Eigen::Vector4d& vec, double r, double theta, double phi) {
    double sigma = Sigma(r, theta);
    double delta = Delta(r, theta);
    double lambda = Lambda(r, theta);

    double p0 = std::sqrt(lambda/(delta*sigma))*vec(0) + 2.0*a*r/std::sqrt(lambda*sigma*delta)*vec(3);
    double p1 = std::sqrt(delta/sigma)*vec(1);
    double p2 = std::sqrt(1.0/sigma)*vec(2);
    double p3 = std::sqrt(sigma/lambda)*1.0/std::sin(theta)*vec(3);

    return Eigen::Vector4d(p0, p1, p2, p3);
}