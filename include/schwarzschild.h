#pragma once

#include <Eigen/Dense>
#include <cmath>

#include "metric.h"

class Schwarzschild : public Metric {
    
private:
    double r_EH = 2.0;

public:
    inline double g_tt(double r, double theta);
    
    inline double g_rr(double r, double theta);

    /*
    inline double g_thth(double r, double theta);

    inline double g_phph(double r, double theta);

    inline double g_tph(double r, double theta);
    */

    void geodesic_eq_rhs(const Vector8d& y, Vector8d& derivatives);

    bool break_integration(const Vector8d& y, bool &outside_celestial_sphere, bool &below_EH);
};