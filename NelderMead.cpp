/*
 * NelderMead.cpp
 *
 *  Created on: May 10, 2011
 *      Author: kyleklein
 */
#include "NelderMead.hpp"
#include <iostream>
#include "string.h"
#include <algorithm>

NelderMead::NelderMead(double *guess, double step, int dimension,
		double (*objFunction)(double *vector, int dimension)) {
	indices = new int[dimension + 1];
	for (int i = 0; i < dimension + 1; i++) {
		indices[i] = i;
	}
	this->simplex = new double[dimension * (dimension + 1)];
	for (int i = 0; i < dimension + 1; i++) {
		for (int j = 0; j < dimension; j++) {
			SIMPLEX(i,j) = guess[j];
			if (i == j + 1)
				SIMPLEX(i,j) += step;
		}
	}
	this->dimension = dimension;
	this->obj_function = objFunction;
	M = new double[dimension];
	obj_function_results = new double[dimension + 1];
	AR = new double[dimension];
	AE = new double[dimension];
	AC = new double[dimension];
	alpha = ALPHA;
	beta = BETA;
	gamma = GAMMA;
	tau = TAU;
}

NelderMead::NelderMead(int dimension,
		double (*objFunction)(double *vector, int dimension)) {
	indices = new int[dimension + 1];
	for (int i = 0; i < dimension + 1; i++) {
		indices[i] = i;
	}
	this->simplex = new double[dimension * (dimension + 1)];
	for (int i = 0; i < dimension + 1; i++) {
		for (int j = 0; j < dimension; j++) {
			SIMPLEX(i,j) = 1.0;
			if (i == j + 1)
				SIMPLEX(i,j) += 1;
		}
	}
	this->dimension = dimension;
	this->obj_function = objFunction;
	M = new double[dimension + 1];
	obj_function_results = new double[dimension + 1];
	AR = new double[dimension];
	AE = new double[dimension];
	AC = new double[dimension];
	alpha = ALPHA;
	beta = BETA;
	gamma = GAMMA;
	tau = TAU;
}

NelderMead::~NelderMead() {
	delete simplex;
	delete M;
	delete obj_function_results;
	delete AR;
	delete AE;
	delete AC;
}

double* NelderMead::solve(int max_iter) {
	//Compute objective function
	for (int i = 0; i < dimension + 1; i++) {
		obj_function_results[i] = obj_function(&SIMPLEX(i,0), dimension);
	}
	sort_simplex(); //Sort the simplex
	double best = obj_function_results[indices[0]];
	int iter = 0;
	double worst = obj_function_results[indices[dimension]];

	while (best > 1e-6 && (max_iter <= 0 || iter < max_iter)) {
		/*if (iter % 250 == 0) {
			std::cout << iter << " " << best << " "
					<< obj_function_results[indices[dimension]] << std::endl;
		}*/
		reflection(); //Compute reflection
		fAR = obj_function(AR, dimension); //Evaluate reflection
		// Case 1: 
		// reflection is better than the best
		if (fAR < best) {
			expansion();
			fAE = obj_function(AE, dimension); //Evaluate expansion

			// Can we accept the expanded point?
			// If expansion is better, use that
			if (fAE < fAR) {
				memmove(&SIMPLEX(dimension, 0), AE, dimension * sizeof(double));
				obj_function_results[indices[dimension]] = fAE;
			} else { //otherwise use reflection
				memmove(&SIMPLEX(dimension, 0), AR, dimension * sizeof(double));
				obj_function_results[indices[dimension]] = fAR;
			}
			//Case 2: [ fAR > best & fAR < the second worst ]
			// --> accept the point
		} else if (fAR <= obj_function_results[indices[dimension - 1]]) {
			memmove(&SIMPLEX(dimension, 0), AR, dimension * sizeof(double));
			obj_function_results[indices[dimension]] = fAR;
			//Case 3
		} else {
			// %inside contraction
			contraction();
			fAC = obj_function(AC, dimension); //Evaluate contraction
			//If Contraction is better, use it
			if (fAC < worst) {
				memmove(&SIMPLEX(dimension, 0), AC, dimension * sizeof(double));
				obj_function_results[indices[dimension]] = fAC;
			} else { //Otherwise, minimize
				// "shrink"
				minimize();
				//re-evaluate for next iteration
				for (int i = 0; i < dimension + 1; i++) {
					obj_function_results[indices[i]] = obj_function(&SIMPLEX(i,0),
							dimension);
				}
			}

		}
		//Resort the indices, as they might be out of order now
		sort_simplex(); //Sort the simplex
		//Find the new best
		best = obj_function_results[indices[0]];
		worst = obj_function_results[indices[dimension]];
		iter++;
	}
	std::cout << iter << " total iterations\n.";
	return &SIMPLEX(0,0);
}

void NelderMead::reflection() {
	for (int i = 0; i < dimension; i++)
		M[i] = 0.0;
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension; j++) {
			M[j] += SIMPLEX(i, j);
			//DIVIDE after, but might need to worry about overflow!
		}
	}
	for (int i = 0; i < dimension; i++) {
		M[i] /= (dimension); //Divide from earlier, then compute
		AR[i] = M[i] + alpha * (M[i] - SIMPLEX(dimension,i));
	}
}
void NelderMead::expansion() {
	for (int i = 0; i < dimension; i++) {
		AE[i] = AR[i] + gamma * (AR[i] - M[i]);
	}
}
void NelderMead::contraction() {
	double *ATilda;
	if (fAR < obj_function_results[indices[dimension]]) {
		ATilda = AR;
	} else {
		ATilda = &SIMPLEX(dimension, 0);
	}
	for (int i = 0; i < dimension; i++) {
		AC[i] = M[i] + beta * (ATilda[i] - M[i]);
	}

}

void NelderMead::minimize() {
	double *ATilda;
	if (fAR < obj_function_results[indices[dimension]]) {
		ATilda = AR;
	} else {
		ATilda = &SIMPLEX(dimension, 0);
	}
	daxpy(&SIMPLEX(dimension,0), tau, &SIMPLEX(0,0), (1.0 - tau), ATilda,
			dimension);
	for (int i = 1; i < dimension; i++) {
		daxpy(&SIMPLEX(i,0), tau, &SIMPLEX(0,0), (1.0 - tau), &SIMPLEX(i,0),
				dimension);
	}

}

// result = scalar1*a + scalar2*b
void NelderMead::daxpy(double *result, double scalar1, double *a,
		double scalar2, double *b, int length) {
	for (int i = 0; i < length; i++) {
		result[i] = scalar1 * a[i] + scalar2 * b[i];
	}
}

//Debugging purposes
void NelderMead::print_simplex() {
	for (int i = 0; i < dimension + 1; i++) {
		for (int j = 0; j < dimension; j++) {
			std::cout << SIMPLEX(i,j) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void NelderMead::sort_simplex() {
	std::sort(indices, indices + dimension + 1,
			IndexSorter(obj_function_results));
}
