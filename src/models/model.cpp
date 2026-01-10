#include "lattice.h"
#include "model.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

// e^-\beta H
template <typename T>
double Model<T>::compute_single_partition_function() {
    T total_energy = compute_total_energy();
    
    double z = exp(- beta * total_energy);

    return z;
};

template <typename T>
double Model<T>::compute_average_energy(std::vector<T> energies) {
    double avg_e = 0;
    for (auto & i: energies) {
        avg_e += i;
    }
    
    return avg_e / energies.size();

};

template <typename T>
double Model<T>::compute_heat_capacity(std::vector<T> energies) {
    double c_v = 0;
    double avg_e = compute_average_energy(energies);

    for (auto & i: energies) {
        i *= i;        
    }

    double avg_e2 = compute_average_energy(energies);

    c_v = avg_e2 - avg_e * avg_e;

    return c_v;
}
