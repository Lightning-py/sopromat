#include "element.hpp"

void element::setFirstDot(dot& d) { this->d1 = d; }

void element::setSecondDot(dot& d) { this->d2 = d; }

void element::setE(double E) { this->E = E; }

void element::setA(double A) { this->A = A; }

dot element::getFirstDot() { return this->d1; }

dot element::getSecondDot() { return this->d2; }

double element::getE() { return this->E; }

double element::getA() { return this->A; }

double element::getLenght() { return this->lenght; }

Matrix2x2 element::getStiffnessMatrix() { return this->stiffness_matrix; }

Matrix4x4 element::getStiffnessMatrixGlobal() {
    return this->stiffness_matrix_global;
}

Matrix2x4 element::getRotationMatrix() { return this->rotation_matrix; }

Matrix4x2 element::getRotationMatrixTransposed() {
    return this->rotation_matrix_transposed;
}

void element::calculate_lenght() {
    this->lenght =
        std::sqrt(std::pow((d1.x - d2.x), 2) + std::pow((d1.y - d2.y), 2));
}

void element::calculate_stiffness_matrix() {
    double coef = (this->E * this->A) / (this->lenght);

    this->stiffness_matrix(0, 0) = coef;
    this->stiffness_matrix(0, 1) = -1 * coef;
    this->stiffness_matrix(1, 0) = -1 * coef;
    this->stiffness_matrix(1, 1) = coef;
}

void element::calculate_rotation_matrix() {
    double x = (d2.x - d1.x);
    double y = (d2.y - d1.y);

    this->rotation_matrix_transposed(0, 0) = x / this->lenght;
    this->rotation_matrix_transposed(1, 0) = y / this->lenght;
    this->rotation_matrix_transposed(2, 1) = x / this->lenght;
    this->rotation_matrix_transposed(3, 1) = y / this->lenght;

    this->rotation_matrix = this->rotation_matrix_transposed.transpose();
}

void element::calculate_stiffness_matrix_global() {
    this->stiffness_matrix_global = this->rotation_matrix_transposed *
                                    this->stiffness_matrix *
                                    this->rotation_matrix;
}

void element::calculate() {
    this->calculate_lenght();
    this->calculate_stiffness_matrix();
    this->calculate_rotation_matrix();
    this->calculate_stiffness_matrix_global();
}

std::ostream& operator<<(std::ostream& os, element& el) {
    os << "First dot:\t" << el.d1 << std::endl;
    os << "Second dot: \t" << el.d2 << std::endl;

    std::cout << "E:\t" << el.E << std::endl;
    std::cout << "A:\t" << el.A << std::endl;

    std::cout << "Lenght:\t" << el.getLenght() << std::endl;

    std::cout << 'Stiffness matrix:' << std::endl;
    std::cout << el.getStiffnessMatrixGlobal() << std::endl;

    std::cout << "Rotation matrix" << std::endl;
    std::cout << el.getRotationMatrix();

    return os;
}

std::istream& operator>>(std::istream& os, element& el) { return os; }