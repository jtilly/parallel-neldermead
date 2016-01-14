/*
 * NelderMead.hpp
 *
 *  Created on: May 10, 2011
 *      Author: kyleklein
 */

#ifndef NELDERMEAD_HPP_
#define NELDERMEAD_HPP_
#define SIMPLEX(i,j) simplex[((indices[(i)])*dimension) + (j) ]
#define ALPHA (1.0)
#define BETA (.5)
#define GAMMA (1.0)
#define TAU (.5)

class NelderMead {
public:
	/**
	 * Given initial guess, a step, the dimension of the simplex,
	 * and a pointer to an objective function.
	 *
	 * simplex: Array of doubles size dimension*(dimension+1)
	 * dimensions: Dimension of each of the dimension+1 vectors.
	 * objFunction: Pointer to the objective function, takes as argument
	 * 	            a vector and its length, should return a double.
	 */
	NelderMead(double *guess, double step, int dimension,
			double (*objFunction)(double *vector, int dimension));

	/**
	 * Same as above except we initialize the simplex to whatever we choose.
	 */
	NelderMead(int dimension,
			double (*objFunction)(double *vector, int dimension));
	/*
	 * Deletes user passed simplex as well as all allocated memory.
	 */
	~NelderMead();
	/**
	 * Find the point which minimizes the objective function, and return
	 * an array of dimension doubles. User is responsible to free that memory.
	 *
	 * Will stop iterating either when the objective function is less than 10^-6,
	 * or after max_iter iterations, whichever happens first. If max_iter < 0,
	 * then it will iterate until objective function is less than 10^-6.
	 */
	double* solve(int max_iter);
	//Set alpha, otherwise assumed to be ALPHA
	void set_alpha(double alpha);
	//Set beta, otherwise assumed BETA
	void set_beta(double beta);
	//Set gamma, otherwise assumed GAMMA
	void set_gamma(double gamma);
	//Set tau, otherwise assumed to be TAU
	void set_tau(double tau);

private:
	void reflection();
	void expansion();
	void contraction();
	void minimize();
	void daxpy(double *result, double scalar1, double *a, double scalar2,
			double *b, int length);
	void print_simplex();
	bool operator()(int i, int j);
	void sort_simplex();
	double *simplex, *M, *AR, *AE, *AC;
	double *obj_function_results;
	double alpha, beta, gamma, tau, fAR, fAE, fAC;
	int *indices;
	int dimension;
	double (*obj_function)(double *vector, int dimension);

};

class IndexSorter {
public:
	IndexSorter(double *arg) :
			obj_function_results(arg) {
	}
	bool operator()(int i, int j) {
		return obj_function_results[i] < obj_function_results[j];
	}
private:
	double *obj_function_results;
};

#endif /* NELDERMEAD_HPP_ */
