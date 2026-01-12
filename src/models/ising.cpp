#include "models/ising.h"
#include "random.h"
#include <iostream>
#include <string>
#include <vector>

using ivec = std::vector<int>;

// H = -J \sum _{ij} s_i s_j - H \sum_i s_i
int IsingModel::compute_spin_neighbours_term(int index) {
    ivec neigh_array = lattice_obj->get_neighbours_array(index);
    int sign_change = lattice_obj->get_lattice_site(index);
    int neigh_sum = 0;
    
    for (auto& i : neigh_array) {
        neigh_sum+= i;
    }

    return neigh_sum * sign_change;

}

int IsingModel::compute_spin_neighbours_term(ivec indices) {
    ivec neigh_array = lattice_obj->get_neighbours_array(indices);
    int sign_change = lattice_obj->get_lattice_site(indices);
    int neigh_sum = 0;
    
    for (auto& i : neigh_array) {
        neigh_sum+= i;
    }

    return neigh_sum * sign_change;

}

int IsingModel::compute_spin_magnetic_term() {
    int magnetic_sum = 0;
    for (auto &i: lattice_obj->get_lattice()) {
        magnetic_sum += i;
    }
    return magnetic_sum;
}

double IsingModel::compute_total_energy() {
    double total_energy;
    int neigh_energy = 0;
    int magnetic_energy = compute_spin_magnetic_term();
    int N = lattice_obj->get_lattice_total_length();

    for (int i = 0; i < N; i++) {
        neigh_energy += compute_spin_neighbours_term(i);
    }
    
    total_energy = -J * 0.5 * neigh_energy - H * magnetic_energy;
    return total_energy;
}

double IsingModel::compute_magnetisation() {
    int magnetisation = compute_spin_magnetic_term();
    return magnetisation / lattice_obj->get_lattice_total_length();
}

void IsingModel::change_spin_randomly(ivec indices) {
    int index = lattice_obj->get_1d_index(indices);
    old_lattice = lattice_obj->get_lattice();
    ivec &lattice = lattice_obj->get_lattice();
    
    lattice[index] *= -1;
}

double IsingModel::compute_energy_diff_flip() {
    double energy_diff=0;
    ivec random_indices = rng::random_int_array(rng::engine, lattice_obj->get_lattice_dim(), 0, lattice_obj->get_lattice_length() - 1);
    int spin = lattice_obj->get_lattice_site(random_indices);
    energy_diff = 2*J * spin * compute_spin_neighbours_term(random_indices) + 2*H*spin;
    change_spin_randomly(random_indices);

    return energy_diff;
}
