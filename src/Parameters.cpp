#include "Parameters.hpp"

double Parameters::getA0() {
    return this->A0;
}

double Parameters::getA1() {
    return this->A1;
}

double Parameters::getCubeSide() {
    return this->a0;
}

double Parameters::getKsi() {
    return this->ksi;
}

double Parameters::getP() {
    return this->p;
}

double Parameters::getQ() {
    return this->q;
}

Parameters &Parameters::withA0(const double &A0) {
    this->A0 = A0;
    return *this;
}

Parameters &Parameters::withA1(const double &A1) {
    this->A1 = A1;
    return *this;
}

Parameters &Parameters::withCubeSide(const double &a0) {
    this->a0 = a0;
    return *this;
}

Parameters &Parameters::withKsi(const double &ksi) {
    this->ksi = ksi;
    return *this;
}

Parameters &Parameters::withP(const double &p) {
    this->p = p;
    return *this;
}

Parameters &Parameters::withQ(const double &q) {
    this->q = q;
    return *this;
}

Parameters &Parameters::operator+(const Parameters &rhs) {
    this->A0 += rhs.A0;
    this->A1 += rhs.A1;
    this->a0 += rhs.a0;
    this->ksi += rhs.ksi;
    this->p += rhs.p;
    this->q += rhs.q;
}

std::ostream &operator<<(std::ostream &out, const Parameters &param) {
    out << "A0 = " << param.A0 << "; "
        << "A1 = " << param.A1 << "; "
        << "ksi = " << param.ksi << "; "
        << "p = " << param.p << "; "
        << "q = " << param.q << "; "
        << "a0 = " << param.a0 << std::endl;
}