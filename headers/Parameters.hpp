#pragma once

#include <ostream>

class Parameters
{
private:
    double A0;
    double A1;
    double p;
    double ksi;
    double q;
    double a0;
public:
    double getA0();
    double getA1();
    double getP();
    double getKsi();
    double getQ();
    double getCubeSide();
    Parameters &withA0(const double &A0);
    Parameters &withA1(const double &A1);
    Parameters &withP(const double &p);
    Parameters &withKsi(const double &ksi);
    Parameters &withQ(const double &q);
    Parameters &withCubeSide(const double &a0);

    Parameters &operator+(const Parameters &rhs);
    friend std::ostream &operator<<(std::ostream &out, const Parameters &param);
};

