/*
 * NelderMead_Driver.cpp
 *
 * Based on the implementations by Kyle Klein and Jeff Borggaard.
 *
 */

#include <cmath>
#include <iostream>
#include <mpi.h>
#include "LeeWiswall.hpp"
#include "ObjFunction.hpp"

int main(int argc, char **argv) {
    int rank = 0, size = 1;
    int prob_size = 300, max_iterations = -1;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (argc <= 1) {
        std::cerr << "Error: incorrect usage ./execname <prob_size>"
        "<max iterations (optional)>\n";
        exit(1);
    }
    if (argc > 2) {
        max_iterations = atoi(argv[2]);
    } else {
        max_iterations = -1;
    }
    
    prob_size = atoi(argv[1]);
    
    double *guess = new double[prob_size];
    for (int i = 0; i < prob_size; i++)
        guess[i] = -1.0;
    
    DistParNelderMead *solver = new DistParNelderMead(guess, 1.0, prob_size, objFunction1,
                                                      rank, size);
    
    double t1, t2;
    t1 = MPI_Wtime();
    double *answer = solver->solve(max_iterations);
    t2 = MPI_Wtime();
    if (rank == 0) {
        std::cout << "Elapsed time during solve: " << t2 - t1 << std::endl;
        for (int i = 0; i < prob_size - 1; i++) {
            std::cout << answer[i] << ", ";
        }
        std::cout << answer[prob_size - 1] << std::endl;
    }
    delete solver;
    MPI_Finalize();
    return 0;
}
