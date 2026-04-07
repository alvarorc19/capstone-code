#include "models/potts.h"
#include "random.h"
#include <iostream>
#include <string>
#include <vector>

using ivec = std::vector<int>;

int PottsModel::delta_function(int a, int b ){
    return (a == b) ? 1 : 0;
}

// H = -J \sum _{ij} delta(s_i, s_j) - \sum H_j \sum_i delta(s_i, s_j)
// See if it would be appropiate to change it to match ising model or something

int PottsModel::compute_spin_neighbours_term(ivec indices) {
    ivec neigh_array = lattice_obj->get_neighbours_array(indices);
    int lattice_value = lattice_obj->get_lattice_site(indices);
    int neigh_sum = 0;

    for (auto &i: neigh_array) {
        neigh_sum += delta_function(lattice_value, i);
    }

    return neigh_sum;
}

int PottsModel::compute_spin_neighbours_term(int index) {
    ivec neigh_array = lattice_obj->get_neighbours_array(index);
    int lattice_value = lattice_obj->get_lattice_site(index);
    int neigh_sum = 0;

    for (auto &i: neigh_array) {
        neigh_sum += delta_function(lattice_value, i);
    }

    return neigh_sum;
}

// // While computing magnetic field the direction is the number of spin
// // Formula:
// // \sum_{sj} H_{sj} \sum_i \delta(sj,si)
// int PottsModel::compute_spin_magnetic_term(int dim) {
//     int magnetic_sum = 0;
//     
//     for (auto &i: lattice_obj->get_lattice()) {
//         magnetic_sum += delta_function(dim, i);
//     }
//
//     return magnetic_sum;
// }

double PottsModel::compute_total_energy() {
    double total_energy;
    int neigh_energy = 0;
    int magnetic_energy = 0;
    int N = lattice_obj->get_lattice_total_length();

    for (int i = 0; i < N; i++) {
        neigh_energy += compute_spin_neighbours_term(i);
    }

    for (int i = 0; i < q; i++) {
        magnetic_energy += vec_H[i] * compute_spin_magnetic_term(i);
    }
    
    // The 1/2 is due to symmetry of delta function
    total_energy = -J * 0.5 * neigh_energy - magnetic_energy;
    return total_energy;
}

////TODO
//Magnetisation here will be the modulus over all the directions
double PottsModel::compute_magnetisation() {
    int magnetisation = 0;
    ivec lattice = lattice_obj->get_lattice();
    for (auto&i:lattice) {
       magnetisation += i; 
    }
    return magnetisation / lattice_obj->get_lattice_total_length();
}

void PottsModel::change_spin_randomly(ivec indices) {
    int index = lattice_obj->get_1d_index(indices);
    int random_number = rng::random_int_number(rng::engine, 0,q);
    ivec &lattice = lattice_obj->get_lattice();
    
    lattice[index] = random_number;
}

double PottsModel::compute_energy_diff_flip() {
    double energy_diff=0;
    ivec random_indices = rng::random_int_array(rng::engine, lattice_obj->get_lattice_dim(), 0, lattice_obj->get_lattice_length() - 1);
    double value = lattice_obj->get_lattice_site(random_indices);
    energy_diff += J * compute_spin_neighbours_term(random_indices) + vec_H[value]*compute_spin_magnetic_term(value);
    change_spin_randomly(random_indices);
    value = lattice_obj->get_lattice_site(random_indices);
    energy_diff -= J * compute_spin_neighbours_term(random_indices) + vec_H[value]*compute_spin_magnetic_term(value);

    return energy_diff;

}
