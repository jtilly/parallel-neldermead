/*
 * ObjFunction.hpp
 *
 *  Created on: May 17, 2011
 *      Author: kyleklein
 */

#ifndef OBJFUNCTION_HPP_
#define OBJFUNCTION_HPP_
#include <algorithm>
#include <math.h>
#include "DistParNelderMead.hpp"
#include <fstream>
#include <iostream>
#include <string.h>


double objFunction1(double *points, int dimension) {
	double sum = 0;
	for (int i = 0; i < dimension; ++i) {
		sum += std::pow(points[i], 2) / dimension;
	}
	return sum;
}

double objFunction2(double *points, int dimension) {
	double sum = 0;
	for (int i = 0; i < dimension; ++i) {
		sum += std::abs(points[i])/ dimension;
	}
	return sum;
}

#endif
