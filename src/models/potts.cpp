#include "models/potts.h"
#include "potts.h"
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
    ivec neigh_array = lattice_obj.get_neighbous_array(indices);
    int lattice_value = lattice_obj.get_lattice_site(indices);
    int neigh_sum = 0;

    for (auto &i: neigh_array) {
        neigh_sum += delta_function(lattice_value, i);
    }

    return neigh_sum;
}

int PottsModel::compute_spin_neighbours_term(int index) {
    ivec neigh_array = lattice_obj.get_neighbous_array(index);
    int lattice_value = lattice_obj.get_lattice_site(index);
    int neigh_sum = 0;

    for (auto &i: neigh_array) {
        neigh_sum += delta_function(lattice_value, i);
    }

    return neigh_sum;
}

// double PottsModel::compute_magnetic_term(int h_j){
//     double j_mag_term = 0;
//     double magnetic_field = vec_H[h_j];
//     int N = lattice_obj.get_lattice_total_length();
//     for (auto & i : lattice_obj.get_lattice()){
//         j_mag_term += delta_function(i, k);
//     }
//
//     return magnetic_field * j_mag_term;
//
// }

int PottsModel::compute_spin_magnetic_term(ivec indices, int h_j) {
    int spin_mag_term;
    int lattice_spin = lattice_obj.get_lattice_site(indices);

    spin_mag_term = delta_function(h_j, lattice_spin);
    return spin_mag_term;
}

int PottsModel::compute_spin_magnetic_term(int index, int h_j) {
    int spin_mag_term;
    int lattice_spin = lattice_obj.get_lattice_site(index);

    spin_mag_term = delta_function(h_j, lattice_spin);
    return spin_mag_term;
}

// To compute magnetic energy have 2 choices:
//  - Fix H to 1 dimension
//  - Compute the magnetic terms like it was a vector, and could do 
//  some lin. alg. stuff with the vectors!!!! Consider how to do this

// TODO fix it to accomodate an H array or a single value
// Have it by default that the H field is in the q direction
int PottsModel::compute_total_energy() {
    int total_energy;
    int neigh_energy = 0;
    int h_j = q;
    int magnetic_energy = 0;
    int N = lattice_obj.get_lattice_total_length();

    for (int i = 0; i < N; i++) {
        neigh_energy += compute_spin_neighbours_term(i);
        magnetic_energy += compute_spin_magnetic_term(i, h_j);
    }
    
    // The 1/2 is due to symmetry of delta function
    total_energy = -J * (1 / 2) * neigh_energy - H * magnetic_energy;
    return total_energy;
}

//TODO change this one
double PottsModel::compute_magnetisation() {
    int magnetisation = 0;
    ivec lattice = lattice_obj.get_lattice();
    for (auto&i:lattice) {
       magnetisation += i; 
    }
    return magnetisation / lattice_obj.get_lattice_total_length();
}

void PottsModel::change_spin_randomly(ivec indices) {
    int index = lattice_obj.get_1d_index(indices);
    int random_number = rng::random_int_concert(rng::engine, 0,q);
    ivec lattice = lattice_obj.get_lattice();
    
    lattice[index] = random_number;
}

void PottsModel::change_spin_randomly(int index) {
    int random_number = rng::random_int_concert(rng::engine, 0,q);
    ivec lattice = lattice_obj.get_lattice();
    
    lattice[index] = random_number;
}
