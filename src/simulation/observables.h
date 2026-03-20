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
    dvec x_magnetisation;
    dvec y_magnetisation;
    double susceptibility;
    double correlation_length;
    dvec correlation_function_array;
    dvec average_cluster_size;
    // renormalised quantities
    dvec rg_x_magnetisation1;
    dvec rg_y_magnetisation1;
    dvec rg_x_magnetisation2;
    dvec rg_y_magnetisation2;
    dvec rg_x_magnetisation3;
    dvec rg_y_magnetisation3;
    dvec rg_energy1;
    dvec rg_energy2;
    dvec rg_energy3;
};


#endif
