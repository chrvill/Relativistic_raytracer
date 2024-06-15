#include <Eigen/Dense>
#include <cmath>
#include "kerr_newman.h"

inline double KerrNewman::rho2(double r, double theta) {
    return r*r + a*a*std::cos(theta)*std::cos(theta);
}

inline double KerrNewman::Delta(double r, double theta) {
    return r*r - 2.0*r + a*a + r_Q2;
}

inline double KerrNewman::g_tt(double r, double theta) {
    double delta = Delta(r, theta);
    double rho2_ = rho2(r, theta);

    return -delta/rho2_ + a*a*std::sin(theta)*std::sin(theta)/rho2_;
}

inline double KerrNewman::g_rr(double r, double theta) {
    double delta = Delta(r, theta);
    double rho2_ = rho2(r, theta);

    return rho2_/delta;
}

inline double KerrNewman::g_thth(double r, double theta) {
    double rho2_ = rho2(r, theta);

    return rho2_;
}

inline double KerrNewman::g_phph(double r, double theta) {
    double delta = Delta(r, theta);
    double rho2_ = rho2(r, theta);

    double sin4 = std::sin(theta)*std::sin(theta)*std::sin(theta)*std::sin(theta);

    return -a*a*sin4*delta/rho2_ + (r*r + a*a)*(r*r + a*a)*std::sin(theta)*std::sin(theta)/rho2_;
}

inline double KerrNewman::g_tph(double r, double theta) {
    double delta = Delta(r, theta);
    double rho2_ = rho2(r, theta);

    return 2.0*a*(std::sin(theta)*std::sin(theta) - (r*r + a*a));
}

