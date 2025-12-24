#include "random.h"
#include "lattice/doubleLattice.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using dvec = std::vector<double>;

DoubleLattice::DoubleLattice(int L, int dim)
    : Lattice(L, dim){
        lattice = generate_lattice(L, dim);
        neighbours_table = calculate_neighbours_table(L, dim);
}

dvec DoubleLattice::generate_lattice(int L , int dim) {
    dvec temp_lattice = rng::random_angle_array(rng::engine,std::pow(L,dim));
    return temp_lattice;
}

double DoubleLattice::get_lattice_site(ivec indices) {
    int index = get_1d_index(indices);
    return lattice[index];
}

double DoubleLattice::get_lattice_site(int index) {
    return lattice[index];
}

dvec DoubleLattice::get_neighbours_array(ivec indices) {
    int index = get_1d_index(indices);
    int neigh_index;
    dvec neigh_array;
    neigh_array.reserve(2*lattice_dim);

    for (int i = 0; i < 2 * lattice_dim; i++) {
        neigh_index = neighbours_table[(2 * lattice_dim) * index + i];
        neigh_array.emplace_back(lattice[neigh_index]);
    }

    return neigh_array;
}

dvec DoubleLattice::get_neighbours_array(int index) {
    int neigh_index;
    dvec neigh_array;
    neigh_array.reserve(2*lattice_dim);

    for (int i = 0; i < 2 * lattice_dim; i++) {
        neigh_index = neighbours_table[(2 * lattice_dim) * index + i];
        neigh_array.emplace_back(lattice[neigh_index]);
    }

    return neigh_array;
}


