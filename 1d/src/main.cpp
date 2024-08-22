#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include "dot.hpp"
#include "element.hpp"
#include "force.hpp"

#define INPUT "input.txt"

#define E (200)
#define A (0.25)

int main() {
    std::unordered_map<uint, dot> number_to_dot;
    std::vector<element> elements;

    // ---------- INPUT

    std::ifstream file;
    file.open(INPUT);

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
        double degrees, mod;

        file >> number >> degrees >> mod;

        forces((number - 1) * 2) = std::cos(degrees * M_PI / 180) * mod;
        forces((number - 1) * 2 + 1) = std::sin(degrees * M_PI / 180) * mod;
        number_to_dot[number].Force = force{degrees, mod};
    }

    file >> block;  //   ELEMENTS
    file >> count;

    for (int i = 0; i < count; ++i) {
        uint first, second;
        file >> first >> second;

        element el(number_to_dot[first], number_to_dot[second], E, A);
        el.calculate();

        elements.push_back(el);
    }

    // ---------    INPUT CHECK

    std::cout << "Got this data (from file " << INPUT << "):" << std::endl
              << std::endl;

    std::cout
        << "----------------------------------------------------------------"
        << std::endl;

    for (int i = 0; i < elements.size(); ++i) {
        std::cout << "Element number: " << i + 1 << std::endl;
        std::cout << elements[i] << std::endl;

        std::cout << "---------------------------------------------------------"
                     "-------"
                  << std::endl
                  << std::endl;
    }

    // ----------   ALGORITHM
    Eigen::MatrixXd stiffness_matrix =
        Eigen::MatrixXd::Zero(2 * dots_count, 2 * dots_count);

    for (int ind = 0; ind < elements.size(); ++ind) {
        uint first = elements[ind].getFirstDot().number - 1;
        uint second = elements[ind].getSecondDot().number - 1;

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                int gi, gj;

                if (i < 2)
                    gi = first;
                else
                    gi = second;

                if (j < 2)
                    gj = first;
                else
                    gj = second;

                stiffness_matrix(gi * 2 + i % 2, gj * 2 + j % 2) +=
                    elements[ind].getStiffnessMatrixGlobal()(i, j);
            }
        }
    }

    std::cout << "Global stiffness matrix before post-processing:" << std::endl;
    std::cout << stiffness_matrix << std::endl << std::endl;

    std::cout << "Forces: " << std::endl;
    std::cout << forces << std::endl << std::endl;

    for (int i = 1; i <= dots_count; ++i) {
        if (number_to_dot[i].x_fixed) {
            for (int j = 0; j < 2 * dots_count; ++j) {
                stiffness_matrix((i - 1) * 2, j) = 0;
                stiffness_matrix(j, (i - 1) * 2) = 0;
            }

            stiffness_matrix((i - 1) * 2, (i - 1) * 2) = 1;
        }

        if (number_to_dot[i].y_fixed) {
            for (int j = 0; j < 2 * dots_count; ++j) {
                stiffness_matrix((i - 1) * 2 + 1, j) = 0;
                stiffness_matrix(j, (i - 1) * 2 + 1) = 0;
            }

            stiffness_matrix((i - 1) * 2 + 1, (i - 1) * 2 + 1) = 1;
        }
    }

    std::cout << "GLobal stiffness matrix after post-processing:" << std::endl;
    std::cout << stiffness_matrix << std::endl << std::endl;

    std::cout << "Determinant of the stiffness matrix:" << std::endl;
    std::cout << stiffness_matrix.determinant() << std::endl << std::endl;

    std::cout << "Movement vector" << std::endl;

    Eigen::VectorXd answer =
        stiffness_matrix.completeOrthogonalDecomposition().solve(forces);

    for (int i = 1; i <= dots_count; ++i) {
        std::cout << i << " x:\t" << answer((i - 1) * 2) << std::endl;
        std::cout << i << " y:\t" << answer((i - 1) * 2 + 1) << std::endl;
    }
}