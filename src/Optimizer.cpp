#include <Optimizer.hpp>
#include <random>
#include <Const.hpp>
#include <iostream>
#include <algorithm>

double getRandomValue() {
    return rand() % 7;
}

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

void configureInDimmer(std::vector<Dot *> &surfData, Dot dot) {
    auto centerSurfDot = std::find(surfData.begin(), surfData.end(), dot);
    (*centerSurfDot)->isDim = true;
}

void configureOnDimmer(std::vector<Dot *> &surfData, Dot dot) {
    surfData.push_back(new Dot(dot));
}

void Optimizer::initBB() {
    for (int i = 0; i < 7; ++i) {
        Parameters params;
        params.withA0(CONST::A0)
            .withA1(0)
            .withQ(getRandomValue())
            .withKsi(getRandomValue())
            .withP(getRandomValue())
            .withCubeSide(getRandomValue());
        simplex.insert({errorBB(params), params});
    }
}

void Optimizer::initAB() {
    for (int i = 0; i < 7; ++i) {
        Parameters params;
        params.withA0(getRandomValue())
            .withA1(getRandomValue())
            .withQ(getRandomValue())
            .withKsi(getRandomValue())
            .withP(getRandomValue())
            .withCubeSide(CONST::a0);
        simplex.insert({errorAB(params), params});
    }
}

void Optimizer::initAA() {
    for (int i = 0; i < 7; ++i) {
        Parameters params;
        params.withA0(getRandomValue())
            .withA1(getRandomValue())
            .withQ(getRandomValue())
            .withKsi(getRandomValue())
            .withP(getRandomValue())
            .withCubeSide(CONST::a0);
        simplex.insert({errorAA(params), params});
    }
}

double Optimizer::errorBB(const Parameters &params) {
    auto gck = init_sharp(params.getCubeSide());

    double eCoh = countEnergy(*gck, params, CONST::dOrig);
    DeformParams deformParams = countParams(*gck, eCoh, params);

    double error = pow(eCoh - CONST::Ec, 2) / 5 +
        pow(deformParams.C11 - CONST::C11, 2) / 5 +
        pow(deformParams.C12 - CONST::C12, 2) / 5 + 
        pow(deformParams.C44 - CONST::C44, 2) / 5 +
        pow(params.getCubeSide() - CONST::a0, 2) / 5;
    delete gck;
    return error;
}

double Optimizer::errorAB(const Parameters &params) {
    double a0 = params.getCubeSide();
    std::vector<Dot*> *surfData = init_sharp(a0);
    (*surfData)[0]->isDim = true;
    double eSurf = countSolEnergy(*surfData, params, CONST::dOrig);
    double eSol = eSurf - CONST::Ec * 108 - CONST::eCohA + CONST::Ec;
    
    double error = pow(eSol - CONST::Esol, 2);

    delete surfData;
    return error;
}

double Optimizer::errorAA(const Parameters &params) {
    double a0 = params.getCubeSide();
    std::vector<Dot*> *surfData = init_sharp(a0);

    double eSurf = countSurfEnergy(*surfData, params, CONST::dOrig);
    configureInDimmer(*surfData, Dot{1.5 * a0, 1 * a0, 2.5 * a0});
    double eAdatom = countSurfEnergy(*surfData, params, CONST::dOrig);
    configureInDimmer(*surfData, Dot{1 * a0, 1.5 * a0, 2.5 * a0});
    double eDimSurf = countSurfEnergy(*surfData, params, CONST::dOrig);

    double EInDim = eDimSurf - eSurf - 2 * (eAdatom - eSurf);
    delete surfData;

    std::vector<Dot*> *surfOnData = init_sharp(a0);
    configureOnDimmer(*surfOnData, Dot{1 * a0, 1 * a0, 3 * a0});
    double eOnAdatom = countSurfEnergy(*surfOnData, params, CONST::dOrig);
    configureOnDimmer(*surfOnData, Dot{1.5 * a0, 1.5 * a0, 3 * a0});
    double eOnDimSurf = countSurfEnergy(*surfOnData, params, CONST::dOrig);

    double EOnDim = eOnDimSurf - eSurf - 2 * (eOnAdatom - eSurf);
    delete surfOnData;

    double error = pow(EInDim - CONST::EinDim, 2) / 2 + pow(EOnDim - CONST::EonDim, 2) / 2;
    return error;
}

void Optimizer::step(double (* error)(const Parameters &params)) {
    Parameters x_h = simplex.rbegin()->second;
    double f_h = simplex.rbegin()->first;
    // Parameters x_g = (++simplex.rbegin())->second;
    double f_g = (++simplex.rbegin())->first;
    Parameters x_l = simplex.begin()->second;
    double f_l = simplex.begin()->first;

    std::cout << "BEST ERROR: " << f_l << std::endl;
    std::cout << "WORST ERROR: " << f_h << std::endl;

    Parameters x_c;

    for(auto it = ++simplex.rbegin(); it != simplex.rend(); ++it) {
        x_c += it->second;
    }
    x_c /= 6;
    Parameters x_r = (x_c*(1+ref_coef) - x_h*ref_coef);
    double f_r = error(x_r);

    if(f_r < f_l) {
        Parameters x_e = (x_c * (1 - stretch_coef) + x_r * stretch_coef);
        double f_e = error(x_e);

        simplex.erase(--simplex.end());
        if(f_e < f_r) {
            simplex.insert({f_e, x_e});
        } else {
            simplex.insert({f_r, x_r});
        }
    } else if(f_r < f_g) {
        simplex.erase(--simplex.end());
        simplex.insert({f_r, x_r});
    } else if(f_r < f_h) {
        simplex.erase(--simplex.end());
        simplex.insert({f_r, x_r});

        Parameters x_s = (x_c * (1 - compr_coef) + x_h * compr_coef);
        double f_s = error(x_s);
        if(f_s < f_h) {
            simplex.erase(--simplex.end());
            simplex.insert({f_s, x_s});
        } else {
            std::multimap<double, Parameters> tmp;
            simplex.clear();

            for(auto& x_i : tmp) {
                Parameters new_vec = x_l + (x_i.second - x_l) * glob_compr_coef;
                simplex.insert({error(new_vec), new_vec});
            }
        }
    } else {

        Parameters x_s = x_c * (1 - compr_coef) + x_h * compr_coef;
        double f_s = error(x_s);
        if(f_s < f_h) {
            simplex.erase(--simplex.end());
            simplex.insert({f_s, x_s});
        } else {
            std::multimap<double, Parameters> tmp(simplex);
            simplex.clear();

            for(auto& x_i : tmp) {
                Parameters new_vec = x_l + (x_i.second - x_l) * glob_compr_coef;
                simplex.insert({error(new_vec), new_vec});
            }
        }
    }
}

void Optimizer::run() {
    initBB();
    while(simplex.begin()->first > 1e-1 || fabs(simplex.rbegin()->first - simplex.begin()->first) < 1e-2) {
        step(&errorBB);
        std::cout << simplex.begin()->second << std::endl;
    }
    std::cout << simplex.begin()->second << std::endl;
}