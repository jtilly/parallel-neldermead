/*
 * test.cpp runs a bunch of tests using the serial implementation of Nelder Mead
 *
 */

#include <cmath>
#include <iostream>
#include "NelderMead.hpp"
#include <stdlib.h>
#include "ObjFunction.hpp"


void test(double *guess, double *truth, int dimension, double (*obj_function)(double *vector, int dimension)) {

    NelderMead *solver = new NelderMead(guess, 1.0, dimension, obj_function);

    double *answer = solver->solve(10000);

    std::cout << "answer" <<  " " << "truth" << std::endl;
    for (int i = 0; i < dimension; i++) {
        std::cout << answer[i] <<  " " << truth[i] <<  " ";
        std::cout << ((std::abs(answer[i] - truth[i]) < 1e-3) ? "pass" : "fail") <<  std::endl;
    }

    delete solver;

}

int main(int argc, char **argv) {

    // rosenbrock
    {
        int dimension = 2;
        std::cout << std::endl << "===================================" << std::endl;
        std::cout << "rosenbrock" << std::endl;

        double *guess = new double[dimension];
        guess[0] = -1.2;
        guess[1] = 1.0;

        double *truth = new double[dimension];
        truth[0] = 1.0;
        truth[1] = 1.0;

        test(guess, truth, dimension, rosenbrock);
        delete[] guess;
        delete[] truth;
    }

    // himmelblau
    {
        int dimension = 2;
        std::cout << std::endl << "===================================" << std::endl;
        std::cout << "himmelblau" << std::endl;

        double *guess = new double[dimension];
        guess[0] = 1.0;
        guess[1] = 1.0;

        double *truth = new double[dimension];
        truth[0] = 3.0;
        truth[1] = 2.0;

        test(guess, truth, dimension, himmelblau);
        delete[] guess;
        delete[] truth;
    }

    // beale
    {
        int dimension = 2;
        std::cout << std::endl << "===================================" << std::endl;
        std::cout << "beale" << std::endl;

        double *guess = new double[dimension];
        guess[0] = 1.0;
        guess[1] = 1.0;

        double *truth = new double[dimension];
        truth[0] = 3.0;
        truth[1] = 0.5;

        test(guess, truth, dimension, beale);
        delete[] guess;
        delete[] truth;
    }

    // bohach2
    {
        int dimension = 2;
        std::cout << std::endl << "===================================" << std::endl;
        std::cout << "bohach2" << std::endl;

        double *guess = new double[dimension];
        guess[0] = 0.6;
        guess[1] = 1.3;

        double *truth = new double[dimension];
        truth[0] = 0.0;
        truth[1] = 0.0;

        test(guess, truth, dimension, bohach2);
        delete[] guess;
        delete[] truth;
    }

    // extended_rosenbrock
    {
        int dimension = 4;
        std::cout << std::endl << "===================================" << std::endl;
        std::cout << "extended_rosenbrock" << std::endl;

        double *guess = new double[dimension];
        guess[0] = 0.0;
        guess[1] = 0.0;
        guess[2] = 0.0;
        guess[3] = 0.0;

        double *truth = new double[dimension];
        truth[0] = 1.0;
        truth[1] = 1.0;
        truth[2] = 1.0;
        truth[3] = 1.0;

        test(guess, truth, dimension, extended_rosenbrock);
        delete[] guess;
        delete[] truth;
    }

    // powell
    {
        int dimension = 4;
        std::cout << std::endl << "===================================" << std::endl;
        std::cout << "powell" << std::endl;

        double *guess = new double[dimension];
        guess[0] = 3.0;
        guess[1] = -1.0;
        guess[2] = 0.0;
        guess[3] = 1.0;

        double *truth = new double[dimension];
        truth[0] = 0.0;
        truth[1] = 0.0;
        truth[2] = 0.0;
        truth[3] = 0.0;

        test(guess, truth, dimension, powell);
        delete[] guess;
        delete[] truth;
    }

    // goldstein_price
    {
        int dimension = 2;
        std::cout << std::endl << "===================================" << std::endl;
        std::cout << "goldstein_price" << std::endl;

        double *guess = new double[dimension];
        guess[0] = -0.5;
        guess[1] = 0.25;

        double *truth = new double[dimension];
        truth[0] = 0.0;
        truth[1] = -1.0;

        test(guess, truth, dimension, goldstein_price);
        delete[] guess;
        delete[] truth;
    }

    // local
    {
        int dimension = 2;
        std::cout << std::endl << "===================================" << std::endl;
        std::cout << "local" << std::endl;

        double *guess = new double[dimension];
        guess[0] = 1.0;
        guess[1] = 1.0;

        double *truth = new double[dimension];
        truth[0] = 0.2858054412;
        truth[1] = 0.2793263206;

        test(guess, truth, dimension, local);
        delete[] guess;
        delete[] truth;
    }

    return 0;
}
