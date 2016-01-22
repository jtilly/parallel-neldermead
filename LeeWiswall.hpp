/*
 * LeeWiswall.hpp
 *
 * Implements MPI based distributed memory parallel NelderMead simplex method.
 *
 * Based on the implementations by Donghoon Lee and Matthew Wiswall,
 * Kyle Klein, and Jeff Borggaard.
 *
 */

#ifndef LEEWISWALL_HPP_
#define LEEWISWALL_HPP_
#define SIMPLEX(i,j) simplex[((indices[(i)])*dimension) + (j) ]
#define RHO (1.0) // RHO > 0
#define XI (2.0)  // XI  > max(RHO, 1)
#define GAM (0.5) // 0 < GAM < 1
#define SIG (0.5) // 0 < SIG < 1
class LeeWiswall {
public:
    /**
     * Given initial guess, a step, the dimension of the simplex,
     * and a pointer to an objective function.
     *
     * simplex: Array of doubles size dimension*(dimension+1)
     * dimensions: Dimension of each of the dimension+1 vectors.
     * obj_function: Pointer to the objective function, takes as argument
     *              a vector and its length, should return a double.
     */
    LeeWiswall(double *guess, double step, int dimension,
                      double (*obj_function)(double *vector, int dimension), int rank, int size);
    
    /**
     * Same as above except we initialize the simplex to whatever we choose.
     */
    LeeWiswall(int dimension,
                      double (*obj_function)(double *vector, int dimension), int rank, int size);
    /*
     * Deletes user passed simplex as well as all allocated memory.
     */
    ~LeeWiswall();
    /**
     * Find the point which minimizes the objective function, and return
     * an array of dimension doubles. User is responsible to free that memory.
     *
     * Will return answer if less than 1e-6, or if max_iterations > 0, then after
     * max_iterations, whichever comes first.
     */
    double* solve(int max_iterations);
    // Set rho, otherwise assumed to be RHO
    void set_rho(double rho);
    // Set xi, otherwise assumed XI
    void set_xi(double xi);
    // Set gam, otherwise assumed GAM
    void set_gam(double gam);
    // Set sig, otherwise assumed to be SIG
    void set_sig(double sig);
    // Set minimimum improvement to do restart after some number of iterations
    void setRestartCriterion(int iterations, double improvement);
    
private:
    void init(double *guess, double step, int dimension,
              double (*obj_function)(double *vector, int dimension), int rank, int size);
    void centroid();
    void reflection();
    void expansion();
    void outsidecontraction();
    void insidecontraction();
    void minimize();
    void broadcast();
    void daxpy(double *result, double scalar1, double *a, double scalar2,
               double *b, int length);
    void print_simplex();
    void sort_simplex();
    void update(double *vector, int index);
    void evaluate_all();
    double *simplex, *M, *AR, *AE, *AC;
    double *obj_function_results;
    double rho, xi, gam, sig, fAR, fAE, fAC, best;
    int *indices;
    int dimension;
    int rank, size, current_point;
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

#endif /* LEEWISWALL_HPP_ */
