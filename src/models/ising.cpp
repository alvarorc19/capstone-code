#include "models/ising.h"
#include "ising.h"
#include <iostream>
#include <string>
#include <vector>

using ivec = std::vector<int>;

//TODO: see how energy is handeled here
// -  Two options:
//    1. Calculate spin stuff in one and then add them all up
//    2. Calculate per site and then just add all the sites
//       - Problem: need to account for the 1/2 in the neighbours term

// H = -J \sum _{ij} s_i s_j - H \sum_i s_i
int IsingModel::compute_spin_magnetic_term(ivec indices) {
    return lattice_obj.get_lattice_site(indices);
    
}

int IsingModel::compute_spin_magnetic_term(int index) {
    return lattice_obj.get_lattice_site(index);
    
}

int IsingModel::compute_spin_neighbours_term(ivec indices) {
    ivec neigh_array = lattice_obj.get_neighbous_array(indices);
    int sign_change = lattice_obj.get_lattice_site(indices);
    int neigh_sum = 0;
    
    for (auto& i : neigh_array) {
        neigh_sum+= i;
    }

    return neigh_sum * sign_change;

}

int IsingModel::compute_spin_neighbours_term(int index) {
    ivec neigh_array = lattice_obj.get_neighbous_array(index);
    int sign_change = lattice_obj.get_lattice_site(index);
    int neigh_sum = 0;
    
    for (auto& i : neigh_array) {
        neigh_sum+= i;
    }

    return neigh_sum * sign_change;

}


int IsingModel::compute_total_energy() {
    int total_energy;
    int neigh_energy = 0;
    int magnetic_energy = 0;
    int N = lattice_obj.get_lattice_total_length();

    for (int i = 0; i < N; i++) {
        neigh_energy += compute_spin_neighbours_term(i);
        magnetic_energy += compute_spin_magnetic_term(i);
    }
    
    total_energy = -J * (1 / 2) * neigh_energy - H * magnetic_energy;
    return total_energy;
}

double IsingModel::compute_magnetisation() {
    int magnetisation = 0;
    ivec lattice = lattice_obj.get_lattice();
    for (auto&i:lattice) {
       magnetisation += i; 
    }
    return magnetisation / lattice_obj.get_lattice_total_length();
}

void IsingModel::change_spin_randomly(ivec indices) {
    int index = lattice_obj.get_1d_index(indices);
    int random_number = rng::random_int_concert(rng::engine, 0,q);
    ivec lattice = lattice_obj.get_lattice();
    
    lattice[index] = random_number;
}

void IsingModel::change_spin_randomly(int index) {
    int random_number = rng::random_int_concert(rng::engine, 0,q);
    ivec lattice = lattice_obj.get_lattice();
    
    lattice[index] = random_number;
}
