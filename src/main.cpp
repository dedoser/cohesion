#include <Dot.hpp>
#include <Const.hpp>
#include <vector>
#include <Params.hpp>
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <stack>
#include <DeformParams.hpp>

std::vector<Dot*> *init_sharp(double a0) {
    std::vector<Dot*> *res = new std::vector<Dot*>();
    for (double i = 0; i < 3; ++i) {
        for (double j = 0; j < 3; ++j) {
            for (double m = 0; m < 3; ++m) {
                res->push_back(new Dot{i * a0, j * a0, m * a0});
                res->push_back(new Dot{(i + 0.5) * a0, (j + 0.5) * a0, m * a0});
                res->push_back(new Dot{i * a0, (j + 0.5) * a0, (m + 0.5) * a0});
                res->push_back(new Dot{(i + 0.5) * a0, j * a0, (m + 0.5) * a0});
            }
        }
    }
    return res;
}

std::vector<Dot *> *init_dim_sharp(double a0) {
    std::vector<Dot*> *res = init_sharp(a0);
    (*res)[0]->isDim = true;
    // (*res)[1]->isDim = true;
    return res;
}

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

double countDimEnergy(const std::vector<Dot*> &data, double a0, const Params *params, const double d[9]) {
    double sumE = 0;
    double cutOff = CONST::CUT_OFF * a0;

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
                        if (dot1.isDim || dot0.isDim) {
                            curER += (params->A1 * (dist - CONST::solr0) + CONST::solA0) * exp(-CONST::solP * (dist / CONST::solr0 - 1));
                            curEB += pow(CONST::solKsi, 2) * exp(-2 * CONST::solQ * (dist / CONST::solr0 - 1));
                        } else {
                            curER += params->A * exp(-params->p * (dist / params->r0 - 1));
                            curEB += pow(params->ksi, 2) * exp(-2 * params->q * (dist / params->r0 - 1));
                        }
                    }
                }
            }
        }
        sumE += (-sqrt(curEB) + curER);
    }
    return sumE;
}

double countSurfDimEnergy(const std::vector<Dot*> &data, double a0, const Params *params, const double d[9]) {
    double sumE = 0;
    double cutOff = CONST::CUT_OFF * a0;

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
                        curER += (CONST::dimA0) * exp(-CONST::dimP * (dist / CONST::dimr0 - 1));
                        curEB += pow(CONST::dimKsi, 2) * exp(-2 * CONST::dimQ * (dist / CONST::dimr0 - 1));
                    } else if (dot1.isDim || dot0.isDim) {
                        curER += (params->A1 * (dist - CONST::solr0) + CONST::solA0) * exp(-CONST::solP * (dist / CONST::solr0 - 1));
                        curEB += pow(CONST::solKsi, 2) * exp(-2 * CONST::solQ * (dist / CONST::solr0 - 1));
                    } else {
                        curER += params->A * exp(-params->p * (dist / params->r0 - 1));
                        curEB += pow(params->ksi, 2) * exp(-2 * params->q * (dist / params->r0 - 1));
                    }
                }
            }
        }
        sumE += (-sqrt(curEB) + curER);
    }
    
    return sumE;
}

DeformParams countParams(double aCur, double ECur, const Params *params) {
    DeformParams deformParams;

    double v0 = pow(aCur, 3), conts_bar = 1.602, alpha = 0.001, alpha2 = pow(alpha, 2);
    double d_b_plus[9] = {1 + alpha, 0, 0, 0, 1 + alpha, 0, 0, 0, 1 + alpha};
    double d_b_minus[9] = {1 - alpha, 0, 0, 0, 1 - alpha, 0, 0, 0, 1 - alpha};
    double d_c11_plus[9] = {1 + alpha, 0, 0, 0, 1 + alpha, 0, 0, 0, 1};
    double d_c11_minus[9] = {1 - alpha, 0, 0, 0, 1 - alpha, 0, 0, 0, 1};
    double d_c12_plus[9] = {1 + alpha, 0, 0, 0, 1 - alpha, 0, 0, 0, 1};
    double d_c12_minus[9] = {1 - alpha, 0, 0, 0, 1 + alpha, 0, 0, 0, 1};
    double d_c44_plus[9] = {1, alpha, 0, alpha, 1, 0, 0, 0, 1 / (1 - alpha2)};
    double d_c44_minus[9] = {1, -alpha, 0, -alpha, 1, 0, 0, 0, 1 / (1 - alpha2)};

    std::vector<Dot *> *data = init_sharp(aCur);

    double E_B_plus = countDimEnergy(*data, aCur, params, d_b_plus) / data->size();
    double E_B_minus = countDimEnergy(*data, aCur, params, d_b_minus) / data->size();
    double E_c11_plus = countDimEnergy(*data, aCur, params, d_c11_plus) / data->size();
    double E_c11_minus = countDimEnergy(*data, aCur, params, d_c11_minus) / data->size();
    double E_c12_plus = countDimEnergy(*data, aCur, params, d_c12_plus) / data->size();
    double E_c12_minus = countDimEnergy(*data, aCur, params, d_c12_minus) / data->size();
    double E_c44_plus = countDimEnergy(*data, aCur, params, d_c44_plus) / data->size();
    double E_c44_minus = countDimEnergy(*data, aCur, params, d_c44_minus) / data->size();

    double d2_E_B = (E_B_plus - 2 * ECur + E_B_minus) / alpha2;
    double d2_E_c11 = (E_c11_plus - 2 * ECur + E_c11_minus) / alpha2;
    double d2_E_c12 = (E_c12_plus - 2 * ECur + E_c12_minus) / alpha2;
    double d2_E_c44 = (E_c44_plus - 2 * ECur + E_c44_minus) / alpha2;

    deformParams.B = 4 * d2_E_B * conts_bar / (9.0 * v0);
    deformParams.C11 = (d2_E_c11 + d2_E_c12) * conts_bar / v0;
    deformParams.C12 = (d2_E_c11 - d2_E_c12) * conts_bar / v0;
    deformParams.C44 = d2_E_c44 * conts_bar / v0;
    delete data;

    return deformParams;
}

void configureInDimmer(std::vector<Dot *> &surfData, Dot dot) {
    auto centerSurfDot = std::find(surfData.begin(), surfData.end(), dot);
    (*centerSurfDot)->isDim = true;
}

void configureOnDimmer(std::vector<Dot *> &surfData, Dot dot) {
    surfData.push_back(new Dot(dot));
}

int main() {
    Params *params = new Params{CONST::a0 / sqrt(2), CONST::A0, CONST::solA1, CONST::Q, CONST::P, CONST::KSI};
    double d0[9] = {1,0,0,0,1,0,0,0,1};
    double aCur = 3.6;
    double step = 0.1;
    double Ecur = 0;
    std::ofstream energy("res/e.txt");
    std::ofstream a("res/a.txt");
    std::stack<double> stack;


    while (fabs(step) > 0.0001) {
        aCur += step;
        std::vector<Dot*> *data = init_sharp(aCur);
        Ecur = countDimEnergy(*data, aCur, params, d0);
        energy << Ecur / data->size() << ' ';
        a << aCur << ' ';
        delete data;
        if (!stack.empty()) {
            double Eprev = stack.top();
            if (Eprev <= Ecur) {
                step = -step / 10;
            }
            stack.push(Ecur);
        } else {
            stack.push(Ecur);
        }
    }
    energy.close();
    a.close();

    Ecur /= 108;

    std::cout << "a0 = " << aCur << "; E0 = " << Ecur << std::endl;

    DeformParams deformParams = countParams(aCur, Ecur, params);

    std::cout << "B = " << deformParams.B
              << "; C11 = " << deformParams.C11
              << "; C12 = " << deformParams.C12
              << "; C44 = " << deformParams.C44
              << std::endl;

    std::vector<Dot*> *data = init_dim_sharp(aCur);
    double eAB = countDimEnergy(*data, aCur, params, d0);
    std::cout << eAB << ' ' << Ecur *  108 << std::endl;
    double eSol = eAB - Ecur * 108 - CONST::eCohA + Ecur;
    std::cout << eSol << std::endl;
    
    delete data;

    std::vector<Dot*> *surfData = init_sharp(aCur);
    double eSurf = countSurfDimEnergy(*surfData, aCur, params, d0);
    std::cout << "Esurf = " << eSurf << std::endl;
    configureInDimmer(*surfData, Dot{1.5 * aCur, 1 * aCur, 2.5 * aCur});
    double eAdatom = countSurfDimEnergy(*surfData, aCur, params, d0);
    std::cout << "E_AdatomSurf = " << eAdatom << std::endl;
    configureInDimmer(*surfData, Dot{1 * aCur, 1.5 * aCur, 2.5 * aCur});
    double eDimSurf = countSurfDimEnergy(*surfData, aCur, params, d0);
    std::cout << "E_DimSurf = " << eDimSurf << std::endl;

    std::cout << "E_in_dim = " << eDimSurf - eSurf - 2 * (eAdatom - eSurf) << std::endl;
    delete surfData;

    std::vector<Dot*> *surfOnData = init_sharp(aCur);
    configureOnDimmer(*surfOnData, Dot{1 * aCur, 1 * aCur, 3 * aCur});
    double eOnAdatom = countSurfDimEnergy(*surfOnData, aCur, params, d0);
    configureOnDimmer(*surfOnData, Dot{1.5 * aCur, 1.5 * aCur, 3 * aCur});
    double eOnDimSurf = countSurfDimEnergy(*surfOnData, aCur, params, d0);
    // for (unsigned i = 0; i < surfOnData->size(); ++i) {
    //     std::cout << *(*surfOnData)[i] << std::endl;
    // }

    std::cout << "E_on_dim = " << eOnDimSurf - eSurf - 2 * (eOnAdatom - eSurf) << std::endl;

    delete params;
    return 0;
}
