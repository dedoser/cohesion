#include "Parameters.hpp"

double Parameters::getA0() const {
    return this->A0;
}

double Parameters::getA1() const {
    return this->A1;
}

double Parameters::getCubeSide() const {
    return this->a0;
}

double Parameters::getKsi() const {
    return this->ksi;
}

double Parameters::getP() const {
    return this->p;
}

double Parameters::getQ() const {
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

Parameters &Parameters::operator+=(const Parameters &rhs) {
    this->A0 += rhs.A0;
    this->A1 += rhs.A1;
    this->a0 += rhs.a0;
    this->ksi += rhs.ksi;
    this->p += rhs.p;
    this->q += rhs.q;
    return *this;
}

Parameters Parameters::operator+(const Parameters &rhs) {
    Parameters parameters = *this;
    parameters += rhs;
    return parameters;
}

Parameters Parameters::operator*(const double &dig) {
    Parameters params = *this;
    params.A0 *= dig;
    params.A1 *= dig;
    params.a0 *= dig;
    params.ksi *= dig;
    params.p *= dig;
    params.q *= dig;
    return params;
}

Parameters Parameters::operator-(const Parameters &rhs) {
    Parameters params = *this;
    params.A0 -= rhs.A0;
    params.A1 -= rhs.A1;
    params.a0 -= rhs.a0;
    params.ksi -= rhs.ksi;
    params.p -= rhs.p;
    params.q -= rhs.q;
    return params;
}

Parameters &Parameters::operator/=(const double &dig) {
    this->A0 /= dig;
    this->A1 /= dig;
    this->a0 /= dig;
    this->ksi /= dig;
    this->p /= dig;
    this->q /= dig;
    return *this;
}

std::ostream &operator<<(std::ostream &out, const Parameters &param) {
    out << "A0 = " << param.A0 << "; "
        << "A1 = " << param.A1 << "; "
        << "ksi = " << param.ksi << "; "
        << "p = " << param.p << "; "
        << "q = " << param.q << "; "
        << "a0 = " << param.a0 << std::endl;
    return out;
}