#include <Dot.hpp>


Dot& Dot::operator-=(const Dot& o) {
    this->x -= o.x;
    this->y -= o.y;
    this->z -= o.z;
    return *this;
}

Dot& Dot::operator*=(double o) {
    this->x *= o;
    this->y *= o;
    this->z *= o;
    return *this;
}

Dot operator-(const Dot& lhs, const Dot& rhs) {
    Dot res(lhs);

    return res -= rhs;
}

Dot operator*(const Dot& lhs, double rhs) {
    Dot res(lhs);

    return res *= rhs;
}

Dot operator*(double lhs, const Dot& rhs) {
    Dot res(rhs);

    return res *= lhs;
}

std::ostream &operator<<(std::ostream &out, const Dot &d) {
    out << "Dot: " << d.x << ' ' << d.y << ' ' << d.z;
    return out;
}

bool operator==(const Dot &lhs, const Dot &right) {
    return lhs.x == right.x &&
           lhs.y == right.y &&
           lhs.z == right.z;
}

bool operator==(const Dot *lhs, const Dot &right) {
    return lhs->x == right.x &&
           lhs->y == right.y &&
           lhs->z == right.z;
}