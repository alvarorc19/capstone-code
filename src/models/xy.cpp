#include "models/xy.h"
#include "random.h"
#include <iostream>
#include <string>
#include <omp.h>
#include <vector>
#include <cmath>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

//H = -J \sum_{ij} cos(θ_i - θ_j) - H_1 \sum_i cos(θ_i) - H_2 \sum_i sin(θ_i)

double XYModel::compute_spin_neighbours_term(int index){
    const ivec& neigh_table = lattice_obj->get_neighbours_table();
    const dvec& lattice = lattice_obj->get_lattice();
    const int stride = lattice_obj->get_lattice_dim() << 1;
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
    const int stride = lattice_dim << 1;
    
    // given neighbours adds them depending on probability
    for (int i = 0; i < stride; i++) { // 2d line
        const int neigh_index = neighbours_table[index * stride + i];

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

// TODO: This thing will only work for 2d at the moment
void XYModel::compute_reduced_lattice(int b, size_t N, size_t L){
    rg_x_spins.clear();
    rg_y_spins.clear();
    const double norm =  1.0 / (b*b);
    rg_x_spins.reserve(N * norm); //2d
    rg_y_spins.reserve(N * norm); // 2d
    const dvec & lattice = lattice_obj->get_lattice();
    int tmp_diag;
    int tmp1;
    int tmp2;
    
    // 2d loop
    for(int i = 0; i < L; i+= b){
        for(int j = 0; j < L; j+= b){
            double x_term = 0;
            double y_term = 0;
            for (int k = 0; k < b; k++){
                tmp_diag = (i+k) * L + (j+k);
                x_term+= std::cos(lattice[tmp_diag]);
                y_term+= std::sin(lattice[tmp_diag]);
                int l = 0;
                while (l < k){
                    tmp1 = (i+k) * L + (j+l);
                    tmp2 = (i+l) * L + (j+k);
                    x_term += std::cos(lattice[tmp1]);
                    x_term += std::cos(lattice[tmp2]);
                    y_term += std::sin(lattice[tmp1]);
                    y_term += std::sin(lattice[tmp2]);
                    l++;
                }
            }
            rg_x_spins.emplace_back(x_term * norm);   
            rg_y_spins.emplace_back(y_term * norm);   
        }
    }

}

void XYModel::compute_reduced_lattice_3d(int b, size_t N, size_t L){
    rg_x_spins.clear();
    rg_y_spins.clear();
    const double norm =  1.0 / (b * b * b);
    rg_x_spins.reserve(N * norm); // 3d
    rg_y_spins.reserve(N * norm); // 3d
    const dvec & lattice = lattice_obj->get_lattice();
    
    // 3d loop
    for(int i = 0; i < L; i+= b){
        for(int j = 0; j < L; j+= b){
            for(int k = 0; k < L; k+=b){
                double x_term = 0;
                double y_term = 0;
                for (int di = 0; di < b; di++) {
                    int row = (i + di) * L * L;
                    for(int dj = 0; dj < b; dj++) {
                        int col = (j + dj) * L;
                        for (int dk = 0; dk < b; dk++){
                            double theta = lattice[row + col + k+dk];
                            x_term += std::cos(theta);
                            y_term += std::sin(theta);
                        }
                    }
                }

                rg_x_spins.emplace_back(x_term * norm);   
                rg_y_spins.emplace_back(y_term * norm);   
            }
        }
    }

}

double XYModel::compute_rg_spin_magnetic_term(int dim){
    double magnetisation = 0;

    if (dim == 0) {
        for (auto &i: rg_x_spins){
            magnetisation += i;
        }
    }

    else if (dim == 1){
        for (auto &i: rg_y_spins){
            magnetisation += i;
        }
    }

    else{
        std::cout << "This only has up to 2 dimensions" << std::endl;
    }

    return magnetisation;
}

// H = -J \sum_<ij> vec(s_i) · vec(s_j) = -J \sum_<ij> s_ix s_jx + s_iy s_jy
double XYModel::compute_rg_energy(int b, size_t N, ivec & neigh_table, int dim){
    double energy = 0;
    const int rg_N = N / ((std::pow(b,dim)));
    const int stride = dim << 1;
    for(int i = 0; i < rg_N; i++) { 
        for (int j = 0; j < stride; j++){
            energy+= rg_x_spins[i] * rg_x_spins[neigh_table[i * stride + j]];
            energy+= rg_y_spins[i] * rg_y_spins[neigh_table[i * stride + j]];
        }
    }
    return -J * (energy / 2);

}

ivec XYModel::calculate_reduced_neighbours_table(int L, int dim, int b) {
    assert(L % b == 0 && "L/b needs to be an exact division to partition it");
    int rg_L = L / b;
    ivec neighbours_array(std::pow(rg_L, dim) * 2 * dim);
    if (dim == 2) {
        for (int i = 0; i < rg_L; i++) {
            for (int j = 0; j < rg_L; j++) {

                int index = i * rg_L + j;

                // i
                int i_up = (i + 1 == rg_L) ? 0 : i + 1;
                int i_down = (i - 1 == -1) ? rg_L - 1 : i - 1;
                // j
                int j_up = (j + 1 == rg_L) ? 0 : j + 1;
                int j_down = (j - 1 == -1) ? rg_L - 1 : j - 1;

                // Up
                neighbours_array[4 * index + 0] = i_up * rg_L + j;
                // Down
                neighbours_array[4 * index + 1] = i_down * rg_L + j;
                // Right
                neighbours_array[4 * index + 2] = i* rg_L + j_up;
                // Left
                neighbours_array[4 * index + 3] = i * rg_L + j_down;
            }
        }

    } else if (dim == 3) {
        for (int i = 0; i < rg_L; i++) {
            for (int j = 0; j < rg_L; j++) {
                for (int k = 0; k < rg_L; k++) {

                    int index = i * rg_L * rg_L + j * rg_L + k;

                    // i
                    int i_up = (i + 1 == rg_L) ? 0 : i + 1;
                    int i_down = (i - 1 == -1) ? rg_L - 1 : i - 1;
                    // j
                    int j_up = (j + 1 == rg_L) ? 0 : j + 1;
                    int j_down = (j - 1 == -1) ? rg_L - 1 : j - 1;
                    // k
                    int k_up = (k + 1 == rg_L) ? 0 : k + 1;
                    int k_down = (k - 1 == -1) ? rg_L - 1 : k - 1;

                    // i
                    neighbours_array[6 * index + 0] = i_up * rg_L * rg_L + j * rg_L + k;
                    neighbours_array[6 * index + 1] = i_down * rg_L * rg_L + j * rg_L + k;
                    // j
                    neighbours_array[6 * index + 2] = i * rg_L * rg_L + j_up * rg_L + k;
                    neighbours_array[6 * index + 3] = i * rg_L * rg_L + j_down * rg_L + k;
                    // k
                    neighbours_array[6 * index + 4] = i * rg_L * rg_L + j * rg_L + k_up;
                    neighbours_array[6 * index + 5] = i * rg_L * rg_L + j * rg_L + k_down;
                }
            }
        }

    } else {
        // TODO Look at obsidian to know how to further implement it
        std::cout << "Further dimensions not implemented yet" << "\n";
        exit(2);
    }
    return neighbours_array;
}




