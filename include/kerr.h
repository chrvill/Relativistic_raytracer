#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <vector>
#include "metric.h"

class Kerr : public Metric {
public:
    double a;
    double r_EH;

    Kerr(double a) : a(a), r_EH(1 + std::sqrt(1 - a * a)) {}

    inline double Sigma(double r, double theta);

    inline double Delta(double r, double theta);

    inline double Lambda(double r, double theta);

    inline double g_tt(double r, double theta);

    inline double g_tph(double r, double theta);

    inline double g_rr(double r, double theta);

    inline double g_thth(double r, double theta);

    inline double g_phph(double r, double theta);

    void geodesic_eq_rhs(const Vector8d& y, Vector8d& derivatives);

    inline bool break_integration(const Vector8d& y, bool &outside_celestial_sphere, bool &below_EH) {
        below_EH = y(1) <= 1.01*r_EH;
        outside_celestial_sphere = y(1) >= 500.0;
        //inside_disk = (std::abs(z) <= 0.1) and (y(1) >= 5.0) and (y(1) <= 10.0);

        return (outside_celestial_sphere or below_EH);
    }
    
    Eigen::Vector3d CartesianToBLVector(const Eigen::Vector3d& vector, double r, double theta, double phi, double a);

    Eigen::Vector3d pos_to_cartesian(double r, double theta, double phi);

    Eigen::Matrix3d transformationMatrix(double r, double theta, double phi);

    Eigen::Vector3d transform_vec_to_cartesian(const Eigen::Vector3d& vec, double r, double theta, double phi);

    Eigen::Vector3d compute_local_cartesian_velocity(const Eigen::Vector4d& u, double r, double theta, double phi);

    Eigen::Vector4d transform_vec_to_global(const Eigen::Vector4d& vec, double r, double theta, double phi);
};
