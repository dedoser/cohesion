#pragma once

#include <ostream>

struct Dot {
    double x, y, z;
    bool isDim = false;
    Dot& operator-= (const Dot& o);
    Dot& operator*= (double o);
    Dot friend operator-(const Dot& lhs, const Dot& rhs);
    Dot friend operator*(const Dot& lhs, double rhs);
    Dot friend operator*(double lhs, const Dot& rhs);
    friend bool operator==(const Dot &lhs, const Dot &right);
    friend bool operator==(const Dot *lhs, const Dot &right);
    friend std::ostream  &operator<<(std::ostream &out, const Dot &d);
};
