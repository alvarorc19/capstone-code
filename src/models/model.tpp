#include "lattice/lattice.h"
#include "models/model.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

// e^-\beta H
template <typename T>
double Model<T>::compute_single_partition_function() {
    T total_energy = compute_total_energy();
    
    double z = exp(- beta * total_energy);

    return z;
}

template <typename T>
double Model<T>::compute_average_energy(std::vector<T> energies) {
    double avg_e = 0;
    for (auto & i: energies) {
        avg_e += i;
    }
    
    return avg_e / energies.size();

}

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

template <typename T>
void Model<T>::write_lattice(HighFive::DataSet* lattice_set, int time) {
    std::vector<T> lattice = lattice_obj->get_lattice();
    size_t N = lattice_obj->get_particle_num();
    std::vector<size_t> new_dims = {static_cast<size_t>(time + 1), static_cast<size_t>(N)};
    lattice_set->resize(new_dims);

    // Generate arrays to select new stuff
    std::vector<size_t> offset = {static_cast<size_t>(time), 0};
    std::vector<size_t> extent = {1, static_cast<size_t>(N)}; 

    // Put lattice in vector
    std::vector<std::vector<T>> data_to_write = {lattice}; 

    // Write the 1xN slice (new row) to the dataset
    lattice_set->select(offset, extent).write(data_to_write);
}
