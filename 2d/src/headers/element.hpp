#ifndef ELEMENT
#define ELEMENT

#include <eigen3/Eigen/Dense>

#include "dot.hpp"
#include "force.hpp"

typedef Eigen::Matrix<double, 3, 3> Matrix3x3;
typedef Eigen::Matrix<double, 8, 8> Matrix8x8;
typedef Eigen::Matrix<double, 3, 8> Matrix3x8;

class element {
   private:
    dot first, second, third, fourth;

    double a, b, c, d, k;

    Matrix8x8 stiffness_matrix = Matrix8x8::Zero();
    Matrix3x3 elasticity_matrix;

   public:
    element(dot& first, dot& second, dot& third, dot& fourth,
            Matrix3x3& elasticity_matrix) {
        this->first = first;
        this->second = second;
        this->third = third;
        this->fourth = fourth;

        this->elasticity_matrix = elasticity_matrix;
    }

    void setFirst(dot& d);
    void setSecond(dot& d);
    void setThird(dot& d);
    void setFourth(dot& d);
    void setStiffnessMatrix(Matrix8x8& matrix);

    dot getFirst();
    dot getSecond();
    dot getThird();
    dot getFourth();
    Matrix8x8 getStiffnessMatrix();

    void calculateCoords();
    void calculateLSM(uint steps = 100);
    Matrix3x8 calculateMatrixB(double x, double y);
};

std::ostream& operator<<(std::ostream& os, element& el);
std::istream& operator>>(std::istream& os, element& el);

#endif  // ELEMENT
