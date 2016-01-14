/*
 * DistParNelderMead.hpp
 *
 * Implements MPI based distributed memory parallel NelderMead simplex method.
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
class DistParNelderMead {
public:
	/**
	 * Given initial guess, a step, the dimension of the simplex,
	 * and a pointer to an objective function.
	 *
	 * simplex: Array of doubles size dimension*(dimension+1)
	 * dimensions: Dimension of each of the dimension+1 vectors.
	 * obj_function: Pointer to the objective function, takes as argument
	 * 	            a vector and its length, should return a double.
	 */
	DistParNelderMead(double *guess, double step, int dimension,
			double (*obj_function)(double *vector, int dimension), int rank, int size,
			int points_per_iter);

	/**
	 * Same as above except we initialize the simplex to whatever we choose.
	 */
	DistParNelderMead(int dimension,
			double (*obj_function)(double *vector, int dimension), int rank, int size,
			int points_per_ter);
	/*
	 * Deletes user passed simplex as well as all allocated memory.
	 */
	~DistParNelderMead();
	/**
	 * Find the point which minimizes the objective function, and return
	 * an array of dimension doubles. User is responsible to free that memory.
	 *
	 * Will return answer if less than 1e-6, or if max_iterations > 0, then after
	 * max_iterations, whichever comes first.
	 */
	double* solve(int max_iterations);
	//Set alpha, otherwise assumed to be ALPHA
	void setAlpha(double alpha);
	//Set beta, otherwise assumed BETA
	void setBeta(double beta);
	//Set gamma, otherwise assumed GAMMA
	void setGamma(double gamma);
	//Set tau, otherwise assumed to be TAU
	void setTau(double tau);
	//Set minimimum improvement to do restart after some number of iterations
	void setRestartCriterion(int iterations, double improvement);

private:
	void init(double *guess, double step, int dimension,
			double (*obj_function)(double *vector, int dimension), int rank, int size,
			int pointsPerIter);
	void centroid();
	void reflection();
	void expansion();
	void contraction();
	void global_best(double *global_best);
	void minimize();
	void daxpy(double *result, double scalar1, double *a, double scalar2,
			double *b, int length);
	void print_simplex();
	void sort_simplex();
	void update(double *vector, int index);
	double *simplex, *M, *AR, *AE, *AC;
	double *obj_function_results;
	double alpha, beta, gamma, tau, fAR, fAE, fAC, best;
	int *indices;
	int dimension, points_on_proc;
	int rank, size, points_per_iter, current_point;
	int updated;
	double (*obj_function)(double *vector, int dimension);

};

class IndexSorter {
public:
	IndexSorter(double *arg) :
			obj_function_results(arg) {};

	bool operator()(int i, int j) {
		return obj_function_results[i] < obj_function_results[j];
	}
private:
	double *obj_function_results;
};

#endif /* NELDERMEAD_HPP_ */
