#ifndef FORCE
#define FORCE

#include <iostream>

typedef struct {
    double x;
    double y;
} force;

std::ostream& operator<<(std::ostream& os, const force& f);
std::istream& operator>>(std::istream& os, force& f);

#endif  // FORCE