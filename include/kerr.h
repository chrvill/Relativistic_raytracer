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

    inline double g_tt(double r, double theta);

    inline double g_tph(double r, double theta);

    inline double g_rr(double r, double theta);

    inline double g_thth(double r, double theta);

    inline double g_phph(double r, double theta);

    Vector8d geodesic_eq_rhs(const Vector8d& y);

    bool break_integration(const Vector8d& y, bool &outside_celestial_sphere, bool &below_EH);

    Eigen::Vector3d CartesianToBLVector(const Eigen::Vector3d& vector, double r, double theta, double phi, double a);

    Eigen::Vector3d pos_to_cartesian(double r, double theta, double phi);

    Eigen::Matrix3d transformationMatrix(double r, double theta, double phi);

    
};
