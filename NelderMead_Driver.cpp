/*
 * NelderMead_Driver.cpp
 *
 *  Created on: May 10, 2011
 *      Author: kyleklein
 */

#include <cmath>
#include <iostream>
#include "NelderMead.hpp"
#include <stdlib.h>
#include "ObjFunction.hpp"


int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Error: incorrect usage ./execname <prob_size>\n";
        exit(1);
    }
    int prob_size = atoi(argv[1]);
    double *guess = new double[prob_size];
    for (int i = 0; i < prob_size; i++)
        guess[i] = -1.0;
    
    NelderMead *solver = new NelderMead(guess, 1.0, prob_size, objFunction3);
    
    double *answer = solver->solve(-1);
    for (int i = 0; i < prob_size - 1; i++) {
        std::cout << answer[i] << ", ";
    }
    std::cout << answer[prob_size - 1] << std::endl;
    delete solver;
    delete guess;
    return 0;
}
