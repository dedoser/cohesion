#include <EnergyFun.hpp>
#include <Const.hpp>
#include <math.h>

Dot countDefformation(const Dot &data, const double d[9]) {
    Dot dot;

    dot.x = data.x * d[0] + data.y * d[1] + data.z * d[2];
    dot.y = data.x * d[3] + data.y * d[4] + data.z * d[5];
    dot.z = data.x * d[6] + data.y * d[7] + data.z * d[8];
    dot.isDim = data.isDim;

    return dot;
}

Dot borderCondition(const Dot &data, const int dx, const int dy, const int dz, const double a0) {
    Dot dot;

    dot.x = data.x + a0 * 3 * dx;
    dot.y = data.y + a0 * 3 * dy;
    dot.z = data.z + a0 * 3 * dz;
    dot.isDim = data.isDim;

    return dot;
}

double countEnergy(const std::vector<Dot*> &data, const Parameters &params, const double d[9]) {
    double sumE = 0;
    double a0 = params.getCubeSide();
    double cutOff = CONST::CUT_OFF * a0;
    double r0 = a0 / sqrt(2);

    #pragma omp parallel for reduction(+: sumE)
    for (unsigned i = 0; i < data.size(); ++i) {
        Dot dot0 = countDefformation(*data[i], d);
        double curEB = 0;
        double curER = 0;
        #pragma omp parallel for reduction(+: curEB, curER)
        for (unsigned j = 0; j < data.size(); ++j) {
            for (int dx = -1; dx < 2; ++dx) {
                for (int dy = -1; dy < 2; ++dy) {
                    for (int dz = -1; dz < 2; ++dz) {
                        if (i == j && dx == 0 && dy == 0 && dz == 0) {
                            continue;
                        }
                        Dot dotBorder = borderCondition(*data[j], dx, dy, dz, a0);
                        Dot dot1 = countDefformation(dotBorder, d);
                        double dist = sqrt(
                            pow(dot1.x - dot0.x, 2) + pow(dot1.y - dot0.y, 2) + pow(dot1.z - dot0.z, 2)
                        );
                        if (dist > cutOff) {
                            continue;
                        }
                        curER += params.getA0() * exp(-params.getP() * (dist / r0 - 1));
                        curEB += pow(params.getKsi(), 2) * exp(-2 * params.getQ() * (dist / r0 - 1));
                    }
                }
            }
        }
        sumE += (-sqrt(curEB) + curER);
    }
    return sumE;
}

DeformParams countParams(std::vector<Dot *> &data, double ECur, const Parameters &params) {
        DeformParams deformParams;

    double v0 = pow(params.getCubeSide(), 3), conts_bar = 1.602, alpha = 0.001, alpha2 = pow(alpha, 2);
    double d_b_plus[9] = {1 + alpha, 0, 0, 0, 1 + alpha, 0, 0, 0, 1 + alpha};
    double d_b_minus[9] = {1 - alpha, 0, 0, 0, 1 - alpha, 0, 0, 0, 1 - alpha};
    double d_c11_plus[9] = {1 + alpha, 0, 0, 0, 1 + alpha, 0, 0, 0, 1};
    double d_c11_minus[9] = {1 - alpha, 0, 0, 0, 1 - alpha, 0, 0, 0, 1};
    double d_c12_plus[9] = {1 + alpha, 0, 0, 0, 1 - alpha, 0, 0, 0, 1};
    double d_c12_minus[9] = {1 - alpha, 0, 0, 0, 1 + alpha, 0, 0, 0, 1};
    double d_c44_plus[9] = {1, alpha, 0, alpha, 1, 0, 0, 0, 1 / (1 - alpha2)};
    double d_c44_minus[9] = {1, -alpha, 0, -alpha, 1, 0, 0, 0, 1 / (1 - alpha2)};

    double E_B_plus = countEnergy(data, params, d_b_plus) / data.size();
    double E_B_minus = countEnergy(data, params, d_b_minus) / data.size();
    double E_c11_plus = countEnergy(data, params, d_c11_plus) / data.size();
    double E_c11_minus = countEnergy(data, params, d_c11_minus) / data.size();
    double E_c12_plus = countEnergy(data, params, d_c12_plus) / data.size();
    double E_c12_minus = countEnergy(data, params, d_c12_minus) / data.size();
    double E_c44_plus = countEnergy(data, params, d_c44_plus) / data.size();
    double E_c44_minus = countEnergy(data, params, d_c44_minus) / data.size();

    double d2_E_B = (E_B_plus - 2 * ECur + E_B_minus) / alpha2;
    double d2_E_c11 = (E_c11_plus - 2 * ECur + E_c11_minus) / alpha2;
    double d2_E_c12 = (E_c12_plus - 2 * ECur + E_c12_minus) / alpha2;
    double d2_E_c44 = (E_c44_plus - 2 * ECur + E_c44_minus) / alpha2;

    deformParams.B = 4 * d2_E_B * conts_bar / (9.0 * v0);
    deformParams.C11 = (d2_E_c11 + d2_E_c12) * conts_bar / v0;
    deformParams.C12 = (d2_E_c11 - d2_E_c12) * conts_bar / v0;
    deformParams.C44 = d2_E_c44 * conts_bar / v0;

    return deformParams;
}

double countSolEnergy(const std::vector<Dot*> &data, const Parameters &params, const double d[9]) {
    double sumE = 0;
    double a0 = params.getCubeSide();
    double cutOff = CONST::CUT_OFF * a0;
    double r0 = a0 / sqrt(2);

    #pragma omp parallel for reduction(+: sumE)
    for (unsigned i = 0; i < data.size(); ++i) {
        Dot dot0 = countDefformation(*data[i], d);
        double curEB = 0;
        double curER = 0;
        #pragma omp parallel for reduction(+: curEB, curER)
        for (unsigned j = 0; j < data.size(); ++j) {
            for (int dx = -1; dx < 2; ++dx) {
                for (int dy = -1; dy < 2; ++dy) {
                    for (int dz = -1; dz < 2; ++dz) {
                        if (i == j && dx == 0 && dy == 0 && dz == 0) {
                            continue;
                        }
                        Dot dotBorder = borderCondition(*data[j], dx, dy, dz, a0);
                        Dot dot1 = countDefformation(dotBorder, d);
                        double dist = sqrt(
                            pow(dot1.x - dot0.x, 2) + pow(dot1.y - dot0.y, 2) + pow(dot1.z - dot0.z, 2)
                        );
                        if (dist > cutOff) {
                            continue;
                        }
                        if (dot0.isDim || dot1.isDim) {
                            curER += (params.getA1() * (dist - r0) + params.getA0()) * exp(-params.getP() * (dist / r0 - 1));
                            curEB += pow(params.getKsi(), 2) * exp(-2 * params.getQ() * (dist / r0 - 1));
                        } else {
                            curER += params.getA0() * exp(-params.getP() * (dist / r0 - 1));
                            curEB += pow(params.getKsi(), 2) * exp(-2 * params.getQ() * (dist / r0 - 1));
                        }
                    }
                }
            }
        }
        sumE += (-sqrt(curEB) + curER);
    }
    return sumE;
}

double countSurfEnergy(const std::vector<Dot*> &data, const Parameters &params, const double d[9]) {
    double sumE = 0;
    double cutOff = CONST::CUT_OFF * params.getCubeSide();
    double a0 = params.getCubeSide();
    double r0 = a0 / sqrt(2);

    #pragma omp parallel for reduction(+: sumE)
    for (unsigned i = 0; i < data.size(); ++i) {
        Dot dot0 = countDefformation(*data[i], d);
        double curEB = 0;
        double curER = 0;
        #pragma omp parallel for reduction(+: curEB, curER)
        for (unsigned j = 0; j < data.size(); ++j) {
            for (int dx = -1; dx < 2; ++dx) {
                for (int dy = -1; dy < 2; ++dy) {
                    if (i == j && dx == 0 && dy == 0) {
                        continue;
                    }
                    Dot dotBorder = borderCondition(*data[j], dx, dy, 0, a0);
                    Dot dot1 = countDefformation(dotBorder, d);
                    double dist = sqrt(
                        pow(dot1.x - dot0.x, 2) + pow(dot1.y - dot0.y, 2) + pow(dot1.z - dot0.z, 2)
                    );
                    if (dist > cutOff) {
                        continue;
                    }
                    if (dot1.isDim && dot0.isDim) {
                        curER += (params.getA1()) * exp(-params.getP() * (dist / r0 - 1));
                        curEB += pow(params.getKsi(), 2) * exp(-2 * params.getQ() * (dist / r0 - 1));
                    } else if (dot1.isDim || dot0.isDim) {
                        curER += (params.getA1() * (dist - r0) + params.getA0()) * exp(-params.getP() * (dist / r0 - 1));
                        curEB += pow(params.getKsi(), 2) * exp(-2 * params.getQ() * (dist / r0 - 1));
                    } else {
                        curER += params.getA0() * exp(-params.getP() * (dist / r0 - 1));
                        curEB += pow(params.getKsi(), 2) * exp(-2 * params.getQ() * (dist / r0 - 1));
                    }
                }
            }
        }
        sumE += (-sqrt(curEB) + curER);
    }
    
    return sumE;
}
