/*
 * NelderMead.cpp
 *
 * Based on the implementations by Kyle Klein and Jeff Borggaard.
 *
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
    rho = RHO;
    xi = XI;
    gam = GAM;
    sig = SIG;
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
            if (i == j + 1) {
                SIMPLEX(i,j) += 1;
            }
        }
    }
    this->dimension = dimension;
    this->obj_function = objFunction;
    M = new double[dimension + 1];
    obj_function_results = new double[dimension + 1];
    AR = new double[dimension];
    AE = new double[dimension];
    AC = new double[dimension];
    rho = RHO;
    xi = XI;
    gam = GAM;
    sig = SIG;
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

        // compute centroid
        centroid();

        // compute reflection and store function value in fAR
        reflection();
        fAR = obj_function(AR, dimension);


        if(best <= fAR && fAR <= obj_function_results[indices[dimension - 1]]) {
            // accept reflection point
            memmove(&SIMPLEX(dimension, 0), AR, dimension * sizeof(double));
            obj_function_results[indices[dimension]] = fAR;
        } else if( fAR < best ) {
            // test for expansion
            expansion();
            fAE = obj_function(AE, dimension);
            if(fAE < fAR) {
                // accept expansion point
                memmove(&SIMPLEX(dimension, 0), AR, dimension * sizeof(double));
                obj_function_results[indices[dimension]] = fAR;
            } else {
                // accept reflection point
                memmove(&SIMPLEX(dimension, 0), AR, dimension * sizeof(double));
                obj_function_results[indices[dimension]] = fAR;
            }
        } else if(obj_function_results[indices[dimension - 1]] <=fAR && fAR < worst) {
            outsidecontraction();
            fAC = obj_function(AC, dimension);

            if(fAC <= fAR) {
                memmove(&SIMPLEX(dimension, 0), AC, dimension * sizeof(double));
                obj_function_results[indices[dimension]] = fAC;
            } else {
                minimize();
                //re-evaluate for next iteration
                for (int i = 0; i < dimension + 1; i++) {
                    obj_function_results[indices[i]] = obj_function(&SIMPLEX(i,0), dimension);
                }
            }
        } else {
            insidecontraction();
            fAC = obj_function(AC, dimension);
            if(fAC < worst) {
                memmove(&SIMPLEX(dimension, 0), AC, dimension * sizeof(double));
                obj_function_results[indices[dimension]] = fAC;
            } else {
                minimize();
                //re-evaluate for next iteration
                for (int i = 0; i < dimension + 1; i++) {
                    obj_function_results[indices[i]] = obj_function(&SIMPLEX(i,0), dimension);
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

void NelderMead::centroid() {
    for (int i = 0; i < dimension; i++) {
        M[i] = 0.0;
    }
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            M[j] += SIMPLEX(i, j);
        }
    }
    for (int i = 0; i < dimension; i++) {
        M[i] /= dimension;
    }
}

void NelderMead::reflection() {
    for (int i = 0; i < dimension; i++) {
        AR[i] = (1 + rho) * M[i] - rho * SIMPLEX(dimension,i);
    }
}
void NelderMead::expansion() {
    for (int i = 0; i < dimension; i++) {
        AE[i] = (1 + rho * xi) * M[i] - rho * xi * SIMPLEX(dimension,i);
    }
}

void NelderMead::insidecontraction() {
    for (int i = 0; i < dimension; i++) {
        AC[i] = (1 - gam) * M[i] + gam * SIMPLEX(dimension,i);
    }
}

void NelderMead::outsidecontraction() {
    for (int i = 0; i < dimension; i++) {
        AC[i] = (1 + rho * gam) * M[i] - rho * gam * SIMPLEX(dimension,i);
    }
}

void NelderMead::minimize() {
    double *ATilda;
    if (fAR < obj_function_results[indices[dimension]]) {
        ATilda = AR;
    } else {
        ATilda = &SIMPLEX(dimension, 0);
    }
    daxpy(&SIMPLEX(dimension,0), sig, &SIMPLEX(0,0), (1.0 - sig), ATilda, dimension);
    for (int i = 1; i < dimension; i++) {
        daxpy(&SIMPLEX(i,0), sig, &SIMPLEX(0,0), (1.0 - sig), &SIMPLEX(i,0), dimension);
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
    std::sort(indices, indices + dimension + 1, IndexSorter(obj_function_results));
}
