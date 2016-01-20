/*
 * NelderMead_Driver.cpp
 *
 * Based on the implementations by Kyle Klein and Jeff Borggaard.
 *
 */

#include <cmath>
#include <iostream>
#include <mpi.h>
#include "DistParNelderMead.hpp"
#include "ObjFunction.hpp"

int main(int argc, char **argv) {
	int rank = 0, size = 1;
	int prob_size = 300, pointsPerIter = 1, max_iterations = -1;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (argc <= 2) {
		std::cerr << "Error: incorrect usage ./execname <prob_size>"
				"<points per iter> <max iterations (optional)>\n";
		exit(1);
	}
	if (argc > 3) {
		max_iterations = atoi(argv[3]);
	} else {
		max_iterations = -1;
	}

	prob_size = atoi(argv[1]);
	pointsPerIter = atoi(argv[2]);
	double *guess = new double[prob_size];
	for (int i = 0; i < prob_size; i++)
		guess[i] = -1.0;
	DistParNelderMead *solver = new DistParNelderMead(guess, 1.0, prob_size, objFunction1,
			rank, size, pointsPerIter);

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
	delete answer;
	MPI_Finalize();
	return 0;
}
