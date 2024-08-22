#include "dot.hpp"

std::ostream& operator<<(std::ostream& os, const dot& d) {
    os << "x:\t" << d.x << std::endl;
    os << "y:\t" << d.y << std::endl;
    os << "x_fixed:\t" << d.x_fixed << std::endl;
    os << "y_fixed:\t" << d.y_fixed << std::endl;
    os << "force:\t" << std::endl << d.Force;

    return os;
}

std::istream& operator>>(std::istream& os, dot& d) {
    os >> d.x >> d.y;

    d.x_fixed = false;
    d.y_fixed = false;

    d.Force = force{0.0, 0.0};

    return os;
}