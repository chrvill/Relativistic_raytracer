#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <vector>
#include "metric.h"

class KerrNewman : public Metric {
public:
    double a;
    double Q;

    double r_Q2;
    double r_EH;

    KerrNewman(double a, double Q): a(a), Q(Q), r_Q2(Q*Q/(4.0*M_PI)), r_EH(1 + std::sqrt(1 - a * a - r_Q2)) {}

    inline double rho2(double r, double theta);
    
    inline double Delta(double r, double theta);

    inline double g_tt(double r, double theta);

    inline double g_tph(double r, double theta);

    inline double g_rr(double r, double theta);

    inline double g_thth(double r, double theta);
    
    inline double g_phph(double r, double theta);

    void geodesic_eq_rhs(const Vector8d& y, Vector8d& derivatives);  

    inline bool break_integration(const Vector8d& y, bool& outside_celestial_sphere, bool& below_EH) {
        below_EH = y(1) <= 1.01*r_EH;
        outside_celestial_sphere = y(1) >= 100.0;

        return (outside_celestial_sphere || below_EH);
    }
};