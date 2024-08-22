#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include "dot.hpp"
#include "element.hpp"
#include "force.hpp"

#define INPUT "input.txt"

// #define INPUT_CHECK
// #define INTEGRATION_CHECK

#define E (200)
#define Nu (0.25)
#define coef (E / (1 - Nu * Nu))

Matrix3x3 global_elasticity_matrix{
    {coef, Nu* coef, 0}, {Nu * coef, coef, 0}, {0, 0, coef * (1 - Nu) / 2}};

int main(int argc, char* argv[]) {
    std::unordered_map<uint, dot> number_to_dot;
    std::vector<element> elements;

    // ---------- INPUT

    std::string INPUT_FILE = INPUT;

    if (argc > 1) {
        std::cout << "Taking file \"" << argv[1] << "\" as input file"
                  << std::endl;
        INPUT_FILE = argv[1];
    }

    std::ifstream file;
    file.open(INPUT_FILE);

    if (!file.is_open()) {
        std::cout << "Error while opening file \"" << INPUT_FILE
                  << "\". Terminating." << std::endl;

        return 1;
    }

    std::string block;  //   DOTS
    file >> block;

    uint dots_count;

    file >> dots_count;
    Eigen::VectorXd forces(2 * dots_count, 1);

    for (int i = 0; i < dots_count; ++i) {
        uint number;
        file >> number;

        dot n;
        n.number = number;
        file >> n;

        number_to_dot[number] = n;
    }

    uint count;

    file >> block;  //   FIXED
    file >> count;

    for (int i = 0; i < count; ++i) {
        uint number;
        char axis;

        file >> number >> axis;

        if (axis == 'x') {
            number_to_dot[number].x_fixed = true;
        } else if (axis == 'y') {
            number_to_dot[number].y_fixed = true;
        } else {
            throw;
        }
    }

    file >> block;  //   FORCE
    file >> count;

    for (int i = 0; i < count; ++i) {
        uint number;

        force f;
        file >> number >> f;

        forces(number - 1) = f.x;
        forces(number - 1 + dots_count) = f.y;

        number_to_dot[number].Force = f;
    }

    file >> block;  //   ELEMENTS
    file >> count;

    for (int i = 0; i < count; ++i) {
        uint first, second, third, fourth;
        file >> first >> second >> third >> fourth;

        element el(number_to_dot[first], number_to_dot[second],
                   number_to_dot[third], number_to_dot[fourth],
                   global_elasticity_matrix);

        el.calculateCoords();
        el.calculateLSM(100);

        elements.push_back(el);
    }

    // ---------    INPUT CHECK

#ifdef INPUT_CHECK
    std::cout << "Got this data (from file " << INPUT << "):" << std::endl
              << std::endl;

    std::cout << "Elements size: " << elements.size() << std::endl;

    for (int i = 0; i < elements.size(); ++i) {
        std::cout << elements[i] << std::endl;
    }
#endif  // INPUT_CHECK

    // ---------    ALGORITHM

    Eigen::MatrixXd stiffness_matrix =
        Eigen::MatrixXd::Zero(2 * dots_count, 2 * dots_count);

    for (int ind = 0; ind < elements.size(); ++ind) {
        element el = elements[ind];

        uint numbers[] = {el.getFirst().number - 1, el.getSecond().number - 1,
                          el.getThird().number - 1, el.getFourth().number - 1};

        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                int indx = 0, indy = 0;

                if (i < 4)  // строки x
                    indx = numbers[i];
                else  // строки y
                    indx = numbers[i - 4] + dots_count;

                if (j < 4)  // столбцы x
                    indy = numbers[j];
                else  // столбцы y
                    indy = numbers[j - 4] + dots_count;

                stiffness_matrix(indx, indy) += el.getStiffnessMatrix()(i, j);
            }
        }
    }

    // ---------- INPUT CHECK

#ifdef INPUT_CHECK
    std::cout << "Global stiffness matrix before post-processing:" << std::endl;
    std::cout << stiffness_matrix << std::endl << std::endl;

    std::cout << "Forces: " << std::endl;
    std::cout << forces << std::endl << std::endl;

#endif  // INPUT_CHECK

    // ---------- ALGORITHM

    for (int i = 1; i <= dots_count; ++i) {
        if (number_to_dot[i].x_fixed) {
            for (int j = 0; j < 2 * dots_count; ++j) {
                stiffness_matrix((i - 1), j) = 0;
                stiffness_matrix(j, (i - 1)) = 0;
            }

            stiffness_matrix((i - 1), (i - 1)) = 1;
        }

        if (number_to_dot[i].y_fixed) {
            for (int j = 0; j < 2 * dots_count; ++j) {
                stiffness_matrix((i - 1) + dots_count, j) = 0;
                stiffness_matrix(j, (i - 1) + dots_count) = 0;
            }

            stiffness_matrix((i - 1) + dots_count, (i - 1) + dots_count) = 1;
        }
    }

    // ---------- INPUT CHECK

#ifdef INPUT_CHECK
    std::cout << "GLobal stiffness matrix after post-processing:" << std::endl;
    std::cout << stiffness_matrix << std::endl << std::endl;

    std::cout << "Determinant of the stiffness matrix:" << std::endl;
    std::cout << stiffness_matrix.determinant() << std::endl << std::endl;

#endif  // INPUT_CHECK

    // ---------- OUTPUT

    std::cout << "Movement vector" << std::endl;

    Eigen::VectorXd answer =
        stiffness_matrix.completeOrthogonalDecomposition().solve(forces);

    for (int i = 1; i <= dots_count; ++i)
        std::cout << i << " x:\t" << answer(i - 1) << std::endl;

    for (int i = 1; i <= dots_count; ++i)
        std::cout << i << " y: \t" << answer(i - 1 + dots_count) << std::endl;
}