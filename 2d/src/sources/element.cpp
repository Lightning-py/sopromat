#include "element.hpp"

inline double getMin(double& first, double& second) {
    if (first <= second) return first;
    return second;
}

inline double getMax(double& first, double& second) {
    if (first >= second) return first;
    return second;
}

inline double getMin4(double& first, double& second, double& third,
                      double& fourth) {
    double temp = first;

    temp = getMin(temp, second);
    temp = getMin(temp, third);
    temp = getMin(temp, fourth);

    return temp;
}

inline double getMax4(double& first, double& second, double& third,
                      double& fourth) {
    double temp = first;

    temp = getMax(temp, second);
    temp = getMax(temp, third);
    temp = getMax(temp, second);

    return temp;
}

void element::setFirst(dot& d) { this->first = d; }

void element::setSecond(dot& d) { this->second = d; }

void element::setThird(dot& d) { this->third = d; }

void element::setFourth(dot& d) { this->fourth = d; }

void element::setStiffnessMatrix(Matrix8x8& matrix) {
    this->stiffness_matrix = matrix;
}

dot element::getFirst() { return this->first; }

dot element::getSecond() { return this->second; }

dot element::getThird() { return this->third; }

dot element::getFourth() { return this->fourth; }

Matrix8x8 element::getStiffnessMatrix() { return this->stiffness_matrix; }

void element::calculateCoords() {
    this->a =
        getMin4(this->first.x, this->second.x, this->third.x, this->fourth.x);
    this->b =
        getMax4(this->first.x, this->second.x, this->third.x, this->fourth.x);

    this->c =
        getMin4(this->first.y, this->second.y, this->third.y, this->fourth.y);
    this->d =
        getMax4(this->first.y, this->second.y, this->third.y, this->fourth.y);

    this->k = 1 / ((a - b) * (c - d));
}

Matrix3x8 element::calculateMatrixB(double x, double y) {
    Matrix3x8 matrixB = Matrix3x8::Zero();

    matrixB(0, 0) = (y - this->d) * k;
    matrixB(0, 1) = -1 * (y - this->d) * k;
    matrixB(0, 2) = (y - this->c) * k;
    matrixB(0, 3) = -1 * (y - this->c) * k;

    matrixB(1, 4) = (x - this->b) * k;
    matrixB(1, 5) = -1 * (x - this->a) * k;
    matrixB(1, 6) = (x - this->a) * k;
    matrixB(1, 7) = -1 * (x - this->b) * k;

    matrixB(2, 0) = (x - this->b) * k;
    matrixB(2, 1) = -1 * (x - this->a) * k;
    matrixB(2, 2) = (x - this->a) * k;
    matrixB(2, 3) = -1 * (x - this->b) * k;
    matrixB(2, 4) = (y - this->d) * k;
    matrixB(2, 5) = -1 * (y - this->d) * k;
    matrixB(2, 6) = (y - this->c) * k;
    matrixB(2, 7) = -1 * (y - this->c) * k;

    return matrixB;
}

void element::calculateLSM(uint steps) {
    double deltaX = (double(this->b - this->a)) / (double(steps));
    double deltaY = (double(this->d - this->c)) / (double(steps));

    double check = 0;

    for (double x = this->a + deltaX / 2; x < this->b; x += deltaX) {
        for (double y = this->c + deltaY / 2; y < this->d; y += deltaY) {
            auto matrixB = this->calculateMatrixB(x, y);

            check += 1;

            this->stiffness_matrix +=
                matrixB.transpose() * this->elasticity_matrix * matrixB;
        }
    }

    check *= (deltaX * deltaY);

#ifdef INTEGRATION_CHECK
    std::cout << "INTEGRATION CHECK " << check << std::endl;
#endif  // INTEGRATION_CHECK

    this->stiffness_matrix *= (deltaX * deltaY);
}

std::ostream& operator<<(std::ostream& os, element& el) {
    os << "First dot:" << std::endl;
    os << el.getFirst() << std::endl;

    os << std::endl << "Second dot:" << std::endl;
    os << el.getSecond() << std::endl;

    os << std::endl << "Third dot:" << std::endl;
    os << el.getThird() << std::endl;

    os << std::endl << "Fourth dot:" << std::endl;
    os << el.getFourth() << std::endl;

    os << "Stiffness matrix:" << std::endl
       << el.getStiffnessMatrix() << std::endl
       << std::endl;

    return os;
}

std::istream& operator>>(std::istream& os, element& el) {
    dot first, second, third, fourth;
    os >> first >> second >> third >> fourth;

    el.setFirst(first);
    el.setSecond(second);
    el.setThird(third);
    el.setFourth(fourth);

    el.calculateCoords();

    return os;
}