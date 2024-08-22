#ifndef DOT
#define DOT

#include <iostream>

#include "force.hpp"

typedef struct {
    uint number;
    double x;
    double y;
    bool x_fixed;
    bool y_fixed;
    force Force;
} dot;

std::ostream& operator<<(std::ostream& os, const dot& d);
std::istream& operator>>(std::istream& os, dot& d);

#endif  // DOT