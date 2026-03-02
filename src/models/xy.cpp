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
    const ivec& neigh_table = lattice_obj->get_neighbours_table();
    const dvec& lattice = lattice_obj->get_lattice();
    const int stride = 2 * lattice_obj->get_lattice_dim();
    double neigh_sum = 0;

    for (int i = 0; i < stride; i++) {
        int neigh_index = neigh_table[index * stride + i];
        neigh_sum += std::cos(lattice[index] - lattice[neigh_index]);
    }

    return neigh_sum;
}
double XYModel::compute_spin_neighbours_term(ivec indices){
    const dvec& neigh_array = lattice_obj->get_neighbours_array(indices);
    const double& lattice_value = lattice_obj->get_lattice_site(indices);
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
    double neigh_energy = 0;
    for (int i = 0; i < lattice_obj->get_particle_num();i++) {
        neigh_energy += compute_spin_neighbours_term(i);
    }

    // // Magnetic term
    // energy+= - vec_H[0] * compute_spin_magnetic_term(0) - vec_H[1] * compute_spin_magnetic_term(1);
    return - J * (neigh_energy / 2);
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
void XYModel::cluster_flip_neighbours(int index, double direction, ivec& cluster_stack, int& spins_flipped, std::vector<uint8_t> & visited, int lattice_dim) {
    const double& index_value = lattice_obj->get_lattice_site(index);
    const ivec& neighbours_table = lattice_obj->get_neighbours_table();
    const dvec& lattice = lattice_obj->get_lattice();
    const double index_dot_prod = 2.0 *  J * beta  * std::cos(direction - index_value);
    
    // given neighbours adds them depending on probability
    for (int i = 0; i < 2 * lattice_dim; i++) {
        const int neigh_index = neighbours_table[index * 2 * lattice_dim + i];

        if (visited[neigh_index]) continue;

        const double neigh_value = lattice[neigh_index];
        const double x =  index_dot_prod * std::cos(direction - neigh_value);

        double probability;

        if ( x <= 0) {
            probability = 1 - exp(x);
        }
        else {
            probability = 0.;
        }

        const double r = rng::random_real_number(rng::engine);
        if (r < probability) {
            flip_spin(neigh_index, direction);
            cluster_stack.emplace_back(neigh_index);
            visited[neigh_index] = 1;
            spins_flipped++;
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
    // old_lattice = lattice;
    lattice[index] = spin;
}

void XYModel::flip_spin(int index, double angle){
    dvec & lattice = lattice_obj->get_lattice();
    // old_lattice = lattice;
    // TODO Hacer esto que te sigue dando error
    double new_angle = M_PI - lattice[index] + 2 * angle;
    if (new_angle < 0){
        lattice[index] = new_angle + 2 * M_PI;
    }
    else if (new_angle > 2 * M_PI){
        lattice[index] = new_angle - 2 * M_PI;
    }
    else {
        lattice[index] = new_angle;
    }
}
