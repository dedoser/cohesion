#pragma once
#include <vector>
#include <Dot.hpp>
#include <Parameters.hpp>
#include <DeformParams.hpp>

DeformParams countParams(std::vector<Dot *> &data, double ECur, const Parameters &params);
double countEnergy(const std::vector<Dot*> &data, const Parameters &params, const double d[9]);
double countSolEnergy(const std::vector<Dot*> &data, const Parameters &params, const double d[9]);
double countSurfEnergy(const std::vector<Dot*> &data, const Parameters &params, const double d[9]);
double countTestEnergy(const std::vector<Dot*> &data, const Parameters &params, const double d[9]);
double countTestSolEnergy(const std::vector<Dot*> &data, const Parameters &params, const double d[9], double a0);
double countTestSurfEnergy(const std::vector<Dot*> &data, const Parameters &params, const double d[9], double a0);