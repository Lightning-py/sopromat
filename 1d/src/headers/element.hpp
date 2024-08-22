#ifndef ELEMENT
#define ELEMENT

#include <eigen3/Eigen/Dense>

#include "dot.hpp"
#include "force.hpp"

typedef Eigen::Matrix<double, 4, 2> Matrix4x2;
typedef Eigen::Matrix<double, 2, 4> Matrix2x4;
typedef Eigen::Matrix2d Matrix2x2;
typedef Eigen::Matrix4d Matrix4x4;

class element {
   private:
    dot d1;
    dot d2;

    double E;
    double A;

    double lenght;

    Matrix2x2 stiffness_matrix;
    Matrix4x4 stiffness_matrix_global;

    Matrix2x4 rotation_matrix = Matrix2x4::Zero();
    Matrix4x2 rotation_matrix_transposed = Matrix4x2::Zero();

   public:
    void setFirstDot(dot& d);
    void setSecondDot(dot& d);
    void setE(double E);
    void setA(double A);

    dot getFirstDot();
    dot getSecondDot();
    double getE();
    double getA();
    double getLenght();
    Matrix2x2 getStiffnessMatrix();
    Matrix4x4 getStiffnessMatrixGlobal();
    Matrix2x4 getRotationMatrix();
    Matrix4x2 getRotationMatrixTransposed();

    void calculate_lenght();
    void calculate_stiffness_matrix();
    void calculate_rotation_matrix();
    void calculate_stiffness_matrix_global();
    void calculate();

    element(dot first, dot second, double E, double A) {
        this->d1 = first;
        this->d2 = second;

        this->E = E;
        this->A = A;

        this->calculate();
    }
};

std::ostream& operator<<(std::ostream& os, element& el);
std::istream& operator>>(std::istream& os, element& el);

#endif  // ELEMENT
