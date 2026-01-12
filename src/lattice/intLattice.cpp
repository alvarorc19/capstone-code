#include "random.h"
#include "lattice/lattice.h"
#include "lattice/intLattice.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using ivec = std::vector<int>;

IntLattice::IntLattice(int L, int dim, int q)
    : Lattice(L, dim), potts_q(q){
        lattice = generate_lattice(L, dim, q);
        neighbours_table = calculate_neighbours_table(L, dim);
}

IntLattice::IntLattice(int L, int dim)
    : Lattice(L, dim){
        lattice = generate_lattice(L, dim);
        neighbours_table = calculate_neighbours_table(L, dim);
}

ivec IntLattice::generate_lattice(int L , int dim) {
    ivec temp_lattice = rng::random_int_array(rng::engine, std::pow(L, dim));

    for (auto &i : temp_lattice) {
        i = (2 * i) - 1;
    }

    return temp_lattice;
}

ivec IntLattice::generate_lattice(int L , int dim, int q) {
    ivec temp_lattice = rng::random_int_array(rng::engine, std::pow(L, dim), 0,q);

    return temp_lattice;
}

int& IntLattice::get_lattice_site(ivec indices) {
    int index = get_1d_index(indices);
    return lattice[index];
}

int& IntLattice::get_lattice_site(int index) {
    return lattice[index];
}

ivec IntLattice::get_neighbours_array(ivec indices) {
    int index = get_1d_index(indices);
    int neigh_index;
    ivec neigh_array;
    neigh_array.reserve(2*lattice_dim);

    for (int i = 0; i < 2 * lattice_dim; i++) {
        neigh_index = neighbours_table[(2 * lattice_dim) * index + i];
        neigh_array.emplace_back(lattice[neigh_index]);
    }

    return neigh_array;
}

ivec IntLattice::get_neighbours_array(int index) {
    int neigh_index;
    ivec neigh_array;
    neigh_array.reserve(2*lattice_dim);

    for (int i = 0; i < 2 * lattice_dim; i++) {
        neigh_index = neighbours_table[(2 * lattice_dim) * index + i];
        neigh_array.emplace_back(lattice[neigh_index]);
    }

    return neigh_array;
}


