/*
 * NelderMead.hpp
 *
 * Based on the implementations by Kyle Klein and Jeff Borggaard.
 *
 */

#ifndef NELDERMEAD_HPP_
#define NELDERMEAD_HPP_
#define SIMPLEX(i,j) simplex[((indices[(i)])*dimension) + (j) ]
#define RHO (1.0) // RHO > 0
#define XI (2.0)    // XI  > max(RHO, 1)
#define GAM (0.5) // 0 < GAM < 1
#define SIG (0.5) // 0 < SIG < 1

class NelderMead {
public:
    /**
     * Given initial guess, a step, the dimension of the simplex,
     * and a pointer to an objective function.
     *
     * simplex: Array of doubles size dimension*(dimension+1)
     * dimensions: Dimension of each of the dimension+1 vectors.
     * objFunction: Pointer to the objective function, takes as argument
     *              a vector and its length, should return a double.
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
    //Set rho, otherwise assumed to be RHO
    void set_rho(double rho);
    //Set xi, otherwise assumed XI
    void set_xi(double xi);
    //Set gam, otherwise assumed GAM
    void set_gam(double gam);
    //Set sig, otherwise assumed to be SIG
    void set_sig(double sig);
    
private:
    void reflection();
    void expansion();
    void insidecontraction();
    void outsidecontraction();
    void minimize();
    void centroid();
    void daxpy(double *result, double scalar1, double *a, double scalar2,
               double *b, int length);
    void print_simplex();
    bool operator()(int i, int j);
    void sort_simplex();
    double *simplex, *M, *AR, *AE, *AC;
    double *obj_function_results;
    double rho, xi, gam, sig, fAR, fAE, fAC;
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
