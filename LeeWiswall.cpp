/*
 * LeeWiswall.cpp
 *
 * Implements MPI based distributed memory parallel NelderMead simplex method.
 *
 * Based on the implementations by Donghoon Lee and Matthew Wiswall,
 * Kyle Klein, and Jeff Borggaard.
 *
 */
#include "LeeWiswall.hpp"
#include <mpi.h>
#include <iostream>
#include "string.h"
#include <algorithm>

DistParNelderMead::DistParNelderMead(double *guess, double step, int dimension,
		double (*obj_function)(double *vector, int dimension), int rank, int size) {
	init(guess, step, dimension, obj_function, rank, size);
}

DistParNelderMead::DistParNelderMead(int dimension,
		double (*obj_function)(double *vector, int dimension), int rank, int size) {
	double *guess = new double[dimension];
	for (int i = 0; i < dimension; i++) {
		guess[i] = 1.0;
	}
	init(guess, 1.0, dimension, obj_function, rank, size);
	delete[] guess;
}

void DistParNelderMead::init(double *guess, double step, int dimension,
		double (*obj_function)(double *vector, int dimension), int rank, int size) {

			indices = new int[(dimension + 1)];
			for (int i = 0; i < (dimension + 1); i++) {
				indices[i] = i;
			}
			this->simplex = new double[dimension * (dimension + 1)];
			for (int i = 0; i < dimension + 1; i++) {
		        for (int j = 0; j < dimension; j++) {
		            SIMPLEX(i,j) = guess[j];
		            if (i == j + 1) {
		                SIMPLEX(i,j) += 1;
		            }
		        }
		    }
			this->dimension = dimension;
			this->obj_function = obj_function;
			this->rank = rank;
			this->size = size;
			M = new double[dimension];
			obj_function_results = new double[(dimension + 1)];
			AR = new double[dimension];
			AE = new double[dimension];
			AC = new double[dimension];
			updated = false;
			rho = RHO;
		    xi = XI;
		    gam = GAM;
		    sig = SIG;
}

DistParNelderMead::~DistParNelderMead() {
	delete[] indices;
	delete[] simplex;
	delete[] M;
	delete[] obj_function_results;
	delete[] AR;
	delete[] AE;
	delete[] AC;
}

double* DistParNelderMead::solve(int max_iterations) {

	//Compute objective function
	for (int i = 0; i < dimension + 1; i++) {
		obj_function_results[i] = obj_function(&SIMPLEX(i, 0), dimension);
	}

	sort_simplex(); //Sort the simplex
	best = obj_function_results[indices[0]];

	int iter = 0;

	while (best > 1e-6 && (max_iterations <= 0 || iter < max_iterations)) {

		updated = 0;

		// which point is this processor replacing?
		current_point = dimension - rank;

		// compute centroid
		centroid();

	    // compute reflection and store function value in fAR
	    reflection();
	    fAR = obj_function(AR, dimension);

	    if(best <= fAR && fAR <= obj_function_results[indices[current_point - 1]]) {
	        // accept reflection point
	        update(AR, current_point);
	        obj_function_results[indices[current_point]] = fAR;
	    } else if(fAR < best) {
	        // test for expansion
	        expansion();
	        fAE = obj_function(AE, dimension);
	        if(fAE < fAR) {
	            // accept expansion point
	            update(AE, current_point);
	            obj_function_results[indices[current_point]] = fAE;
	        } else {
	            // eventual accept reflection point
	            update(AR, current_point);
	            obj_function_results[indices[current_point]] = fAR;
	        }
	    } else if(obj_function_results[indices[current_point - 1]] <=fAR && fAR < obj_function_results[indices[current_point]]) {
	        // do outside contraction
	        outsidecontraction();
	        fAC = obj_function(AC, dimension);
	        if(fAC <= fAR) {
	            // accept outside contraction point
	            update(AC, current_point);
	            obj_function_results[indices[current_point]] = fAC;
	        } else {
	            // shrink
	            memmove(&SIMPLEX(current_point, 0), AR, dimension * sizeof(double));
	            obj_function_results[indices[current_point]] = fAR;
	        }
	    } else {
	        // do inside contraction
	        insidecontraction();
	        fAC = obj_function(AC, dimension);
	        if(fAC < obj_function_results[indices[current_point]]) {
	            // accept inside contraction point
	            update(AC, current_point);
	            obj_function_results[indices[current_point]] = fAC;
	        } else {
	            // shrink
	            memmove(&SIMPLEX(current_point, 0), AR, dimension * sizeof(double));
	            obj_function_results[indices[current_point]] = fAR;
	        }
	    }

		int global_updated = 0;
		MPI_Allreduce(&updated, &global_updated, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if (!global_updated) { //not one processor had an update, minimize
			minimize();
			//Re-eval all of the points
			for (int i = 0; i < dimension + 1; i++) {
				obj_function_results[indices[i]] = obj_function(&SIMPLEX(i, 0), dimension);
			}
		} else {
			broadcast();
		}
		sort_simplex(); //Sort the simplex
		//Find the new best
		best = obj_function_results[indices[0]];

		/*if (iter * size % 500 == 0 && rank == 0) {
				std::cout << iter << " " " " << best << std::endl;
			}*/
		iter++;

	}
	if (rank == 0) {
		std::cout << "Total Iterations: " << iter << std::endl;
	}

	return &SIMPLEX(0,0);
}

void DistParNelderMead::update(double *vector, int index) {
	if (!updated) { //only need to check if not already updated
		for (int i = 0; i < dimension; i++) {
			if (vector[i] != SIMPLEX(index, i)) {
				updated = 1;
				break;
			}
		}
	}
	if (updated) { //might be a new vector, copy it in
		memmove(&SIMPLEX(index, 0), vector, dimension * sizeof(double));
	}
}

void DistParNelderMead::centroid() {
	for (int i = 0; i < dimension; i++) {
		M[i] = 0.0;
	}

	for (int i = 0; i < dimension + 1 - size; i++) {
		for (int j = 0; j < dimension; j++) {
			M[j] += SIMPLEX(i, j);
			//Divide after. Possible overflow for large obj function values!
		}
	}
	for (int i = 0; i < dimension; i++) {
		M[i] /= (dimension + 1 - size); //Divide from earlier, then compute
	}
}

void DistParNelderMead::reflection() {
    for (int i = 0; i < dimension; i++) {
        AR[i] = (1 + rho) * M[i] - rho * SIMPLEX(current_point,i);
    }
}

void DistParNelderMead::expansion() {
    for (int i = 0; i < dimension; i++) {
        AE[i] = (1 + rho * xi) * M[i] - rho * xi * SIMPLEX(current_point,i);
    }
}

void DistParNelderMead::insidecontraction() {
    for (int i = 0; i < dimension; i++) {
        AC[i] = (1 - gam) * M[i] + gam * SIMPLEX(current_point,i);
    }
}

void DistParNelderMead::outsidecontraction() {
    for (int i = 0; i < dimension; i++) {
        AC[i] = (1 + rho * gam) * M[i] - rho * gam * SIMPLEX(current_point,i);
    }
}

void DistParNelderMead::broadcast() {

	double *border_simplex = new double[dimension * size];
	double *border_fval = new double[size];

	MPI_Allgather(&SIMPLEX(current_point, 0), dimension, MPI_DOUBLE, border_simplex, dimension, MPI_DOUBLE,  MPI_COMM_WORLD);
	MPI_Allgather(&(obj_function_results[indices[current_point]]), 1, MPI_DOUBLE, border_fval, 1, MPI_DOUBLE,  MPI_COMM_WORLD);

	// use border_simplex to assemble new simplex
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < dimension; j++) {
			SIMPLEX(dimension + 1 - size + i, j) = border_simplex[i * dimension + j];
		}
	}

	// use border_fval to assemble new fval
	for (int i = 0; i < size; i++) {
		obj_function_results[indices[dimension + 1 - size + i]] = border_fval[i];
	}

}

void DistParNelderMead::minimize() {
	for (int i = 1; i < dimension + 1; i++) {
        daxpy(&SIMPLEX(i,0), sig, &SIMPLEX(0,0), (1.0 - sig), &SIMPLEX(i,0), dimension);
    }
}

// result = scalar1*a + scalar2*b
void DistParNelderMead::daxpy(double *result, double scalar1, double *a,
		double scalar2, double *b, int length) {
	for (int i = 0; i < length; i++) {
		result[i] = scalar1 * a[i] + scalar2 * b[i];
	}
}

//Debugging purposes
void DistParNelderMead::print_simplex() {
	for (int i = 0; i < dimension + 1; i++) {
		for (int j = 0; j < dimension; j++) {
			std::cout << SIMPLEX(i, j) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void DistParNelderMead::sort_simplex() {
	std::sort(indices, indices + dimension + 1, IndexSorter(obj_function_results));
}
