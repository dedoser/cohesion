#pragma once

#include <map>
#include <Parameters.hpp>
#include <EnergyFun.hpp>

class Optimizer
{
private:
    std::multimap<double, Parameters> simplex;
    Parameters BBParams;
    Parameters ABParams;
    Parameters AAParams;
    double ref_coef = 1.0;
    double compr_coef = 0.5;
    double stretch_coef = 2.0;
    double glob_compr_coef = 0.5;
public:
    static double errorBB(const Parameters &params);
    static double errorAB(const Parameters &params);
    static double errorAA(const Parameters &params);
    void initBB();
    void initAB();
    void initAA();
    void step(double (* error)(const Parameters &params));
    void run();
};
