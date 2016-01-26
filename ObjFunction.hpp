/*
 * ObjFunction.hpp
 *
 * Based on the implementations by Kyle Klein and Jeff Borggaard.
 *
 */

#ifndef OBJFUNCTION_HPP_
#define OBJFUNCTION_HPP_
#include <algorithm>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>


// Computes the Himmelblau function.
//
// This function has 4 global minima:
//
//   X* = (  3,        2       ), F(X*) = 0.
//   X* = (  3.58439, -1.84813 ), F(X*) = 0.
//   X* = ( -3.77934, -3.28317 ), F(X*) = 0.
//   X* = ( -2.80512,  3.13134 ), F(X*) = 0.
//
// Suggested starting points are
//
//   (+1,+1),
//   (-1,+1),
//   (-1,-1),
//   (+1,-1),
double himmelblau(double *points, int dimension) {
    
    double sum = 0;
    sum += std::pow(std::pow(points[0], 2) + points[1] - 11, 2);
    sum += std::pow(points[0] + std::pow(points[1], 2) - 7, 2);
    
    return sum;
}

// Computes the Rosenbrock function.
//
// There is a global minimum at X* = (1,1), F(X*) = 0.
//
// The starting point X = [ -1.2, 1.0 ] is suggested.
//
// The contours are sharply twisted.
double rosenbrock(double *points, int dimension) {
    
    double sum = 0;
    sum += std::pow(1.0 - points[0], 2);
    sum += 100.0 * std::pow(points[1] - points[0] * points[0], 2);
    
    return sum;
}

// Computes the Beale function.
//
// This function has a global minimizer:
//
//   X = ( 3.0, 0.5 )
//
// for which
//
//   F(X) = 0.
//
// For a relatively easy computation, start the search at
//
//   X = ( 1.0, 1.0 )
//
// A harder computation starts at
//
//   X = ( 1.0, 4.0 )
double beale(double *points, int dimension) {
    
    double sum = 0;
    sum += std::pow(1.5   - points[0] * (1.0 - points[1]), 2);
    sum += std::pow(2.25  - points[0] * (1.0 - std::pow(points[1], 2)), 2);
    sum += std::pow(2.625 - points[0] * (1.0 - std::pow(points[1], 3)), 2);
    
    return sum;
    
}

// Evaluates the Bohachevsky function #2.
//
// The minimizer is
//
//   X* = [ 0.0, 0.0 ]
//   F(X*) = 0.0
//
// Suggested starting point:
//
//   X init = [ 0.6, 1.3 ];
double bohach2(double *points, int dimension) {
    
    double sum = 0;
    const double pi = 3.14159265358979323846;
    
    sum += points[0] * points[0];
    sum +=  2.0 * points[1] * points[1];
    sum -=  0.3 * cos(3.0 * pi * points[0]) * cos(4.0 * pi * points[1]);
    sum += + 0.3;
    
    return sum;
    
}

//  Computes the extended Rosenbrock function.
//
// The number of dimensions is arbitrary, except that it must be even.
//
// There is a global minimum at X* = (1,1,...), F(X*) = 0.
//
// The contours are sharply twisted.
double extended_rosenbrock(double *points, int dimension) {
    
    double sum = 0;
    double *r = new double[dimension];
    
    if(dimension % 2 > 0) {
        std::cerr << "Dimension must be an even number for extended Rosenbrock function.";
        exit(1);
    }
    
    for (int i = 0; i < dimension; i += 2) {
        r[i] = 1.0 - points[i];
        r[i + 1] = 10.0 * (points[i + 1] - std::pow(points[i], 2));
    }
    
    for (int i = 0; i < dimension; i++) {
        sum += std::pow(r[i], 2);
    }
    
    delete[] r;
    
    return sum;
    
}

// Computes the Powell singular quartic function.
//
// This function has a global minimizer:
//
//   X* = ( 0.0, 0.0, 0.0, 0.0 ), F(X*) = 0.
//
// Start the search at
//
//   X = ( 3.0, -1.0, 0.0, 1.0 )
double powell(double *points, int dimension) {
    
    double f1 = points[0] + 10.0 * points[1];
    double f2 = points[2] - points[3];
    double f3 = points[1] - 2.0 * points[2];
    double f4 = points[0] - points[3];
    
    return f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4;
    
}

//  Evaluates the Goldstein-Price polynomial.
//
//  The minimizer is
//
//  X* = [ 0.0, -1.0 ]
//  F(X*) = 3.0
//
//  Suggested starting point:
//
//  X init = [ -0.5, 0.25 ] (easy convergence)
//  X init = [ -4.0, 5.00 ] (harder convergence)
double goldstein_price(double *points, int dimension) {
    
    double a = points[0] + points[1] + 1.0;
    double b = 19.0 - 14.0 * points[0] + 3.0 * points[0] * points[0] - 14.0 * points[1] + 6.0 * points[0] * points[1] + 3.0 * points[1] * points[1];
    double c = 2.0 * points[0] - 3.0 * points[1];
    double d = 18.0 - 32.0 * points[0] + 12.0 * points[0] * points[0] + 48.0 * points[1] - 36.0 * points[0] * points[1] + 27.0 * points[1] * points[1];
    double f = ( 1.0 + a * a * b ) * ( 30.0 + c * c * d );
    
    return f;
}


//  Computes the local function.
//
//  This function has a local minimizer:
//
//    X* = ( 0.2858054412..., 0.2793263206...), F(X*) = 5.9225...
//
//  and a global minimizer:
//
//    X* = ( -21.02667179..., -36.75997872...), F(X*) = 0.
//
//  Suggested starting point for local minimizer:
//
//    X = ( 1, 1 ), F(X) = 3.33 * 10^6.
//
//  Suggested starting point for global minimizer:
//
//    X = ( -15, -35), F(X) = 1.49 * 10^8.
double local(double *points, int dimension) {
    
    double sum = 0;
    
    sum += std::pow(std::pow(points[0], 2) + 12 * points[1] - 1, 2);
    sum += std::pow(49 * std::pow(points[0], 2) + 49 * std::pow(points[1], 2) + 84 * points[0] + 2324 * points[1] - 681, 2);
    
    return sum;
}

//  Computes the sum of squares function
double objFunction1(double *points, int dimension) {
    
    double sum = 0;
    
    for (int i = 0; i < dimension; ++i) {
        sum += std::pow(points[i] - 2.0, 2) / dimension;
    }
    
    return sum;
}

//  Computes the sum of absolute values function
double objFunction2(double *points, int dimension) {
    
    double sum = 0;
    
    for (int i = 0; i < dimension; ++i) {
        sum += std::abs(points[i] - 2.0)/ dimension;
    }
    
    return sum;
}

//  Computes the sum of squares function
double objFunction3(double *points, int dimension) {
    
    double sum = 0;
    
    for (int i = 0; i < dimension; ++i) {
        sum += std::pow(points[i] - (double)i/2.0, 2) / dimension;
    }
    
    return sum;
}

#endif
