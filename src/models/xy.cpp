#include "models/xy.h"
#include "random.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

//H = -J \sum_{ij} cos(θ_i - θ_j) - H_1 \sum_i cos(θ_i) - H_2 \sum_i sin(θ_i)


double XYModel::compute_spin_neighbours_term(int index){
    dvec neigh_array = lattice_obj->get_neighbours_array(index);
    double lattice_value = lattice_obj->get_lattice_site(index);
    double neigh_sum = 0;

    for (auto &i: neigh_array) {
        neigh_sum += std::cos(lattice_value - i);
    }

    return neigh_sum;
}
double XYModel::compute_spin_neighbours_term(ivec indices){
    dvec neigh_array = lattice_obj->get_neighbours_array(indices);
    double lattice_value = lattice_obj->get_lattice_site(indices);
    double neigh_sum = 0;

    for (auto &i: neigh_array) {
        neigh_sum += std::cos(lattice_value - i);
    }

    return neigh_sum;
}
double XYModel::compute_spin_magnetic_term(int dim){
    double spin_mag_term = 0;
    if (dim == 0) {
       for (auto &i: lattice_obj->get_lattice()){
           spin_mag_term += std::cos(i);
       }
    }
    else if (dim == 1) {
       for (auto &i: lattice_obj->get_lattice()){
           spin_mag_term += std::sin(i);
       }
    }
    else {
        return 0;
    }
    return spin_mag_term;
}
double XYModel::compute_total_energy(){
    double energy = 0;
    // Neighbour term
    double neigh_energy = 0;
    for (int i = 0; i< lattice_obj->get_particle_num();i++) {
        neigh_energy += compute_spin_neighbours_term(i);
    }

    energy +=- J * (neigh_energy / 2);

    // // Magnetic term
    // energy+= - vec_H[0] * compute_spin_magnetic_term(0) - vec_H[1] * compute_spin_magnetic_term(1);
    return energy;
}
//Magnetisation here will be taken as the modulus
// m = sqrt(M_x^2 + M_y^2)
double XYModel::compute_magnetisation(){
    int N = lattice_obj->get_particle_num();
    double M_x = compute_spin_magnetic_term(0) / N;
    double M_y = compute_spin_magnetic_term(1) / N;
    double m = std::sqrt(M_x * M_x + M_y * M_y);
    return m;
    }
// For θi->θ'i, difference is given by
// ΔE = H'-H = -vecH\cdot(vec_s'i-vec_si) - J (sum_j cos(θ'i-θj) - cos(θi - θj))
double XYModel::compute_energy_diff_flip(){
    double energy_diff=0;
    ivec random_indices = rng::random_int_array(rng::engine, lattice_obj->get_lattice_dim(), 0, lattice_obj->get_lattice_length() - 1);
    double angle = lattice_obj->get_lattice_site(random_indices);
    energy_diff += J * compute_spin_neighbours_term(random_indices);

    change_spin_randomly(random_indices);
    angle = lattice_obj->get_lattice_site(random_indices);
    energy_diff -= J * compute_spin_neighbours_term(random_indices);

    return energy_diff;
}

/**
* @brief Adds a spin to the cluster
* 
* @details $$
* P (s_i, s_j) = 1- \exp(\min\{0, \, 2\beta \cos(\phi - \theta_i) \cos(\phi - \theta_j)\})
* $$
* For $\mathbf{r} = (\cos\phi, \sin\phi)$ and $\mathbf{s}_i = (\cos\theta_i, \sin\theta_i)$.
*/
void XYModel::cluster_flip_neighbours(int index, double direction, int& new_spins_flipped, double angle_flip) {
    ivec neigh_indices = lattice_obj->get_neighbours_indices(index);

    // Given neighbours adds them depending on probability
    for (auto & neigh_index: neigh_indices) {
        double index_value = lattice_obj->get_lattice_site(index);
        double neigh_value = lattice_obj->get_lattice_site(neigh_index);
        double probability = 1 - exp(std::min(0., 2 * beta * J * std::cos(direction - index_value) * std::cos(direction - neigh_value)));
        double r = rng::random_real_number(rng::engine);
        if (r < probability) {
            change_spin(neigh_index, angle_flip);
            new_spins_flipped++;
        }
    }
}

void XYModel::change_spin_randomly(ivec indices) {
    int index = lattice_obj->get_1d_index(indices);
    dvec & lattice = lattice_obj->get_lattice();
    old_lattice = lattice;
    lattice[index] = rng::random_angle(rng::engine);
}

void XYModel::change_spin(int index, double spin){
    dvec & lattice = lattice_obj->get_lattice();
    old_lattice = lattice;
    lattice[index] = spin;
}
