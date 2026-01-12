#ifndef OBSERVABLES_H
#define OBSERVABLES_H
#include <iostream>
#include <string>
#include <vector>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

//*********************
// magnetisation
// energy
// susceptibility
// specific heat
// correlation length
// correlation function


struct Observables {
    dvec energy_array;
    dvec magnetisation_array;
    double susceptibility;
    double correlation_length;
    dvec correlation_function_array;
};


#endif
