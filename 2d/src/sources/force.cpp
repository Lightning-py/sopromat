#include "force.hpp"

std::ostream& operator<<(std::ostream& os, const force& f) {
    os << "x:\t" << f.x << "\ny:t" << f.y << std::endl;

    return os;
}

std::istream& operator>>(std::istream& os, force& f) {
    os >> f.x >> f.y;

    return os;
}