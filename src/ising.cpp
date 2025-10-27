#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using ivec = std::vector<int>;
ivec initialise_lattice(int L);
double compute_energy_at_site(int x, int y, ivec lattice);
double compute_total_energy(ivec lattice, double J, double H);
void metropolis_algorithm(ivec &lattice, double beta, double J, double H);
void writer(int iteration, ivec lattice, std::ofstream &file);

int main(int argc, char *argv[]) {
    std::ofstream ofs("results.csv");

    int L = 20;
    double J = -30;
    double H = -200;
    ivec lattice;
    lattice = initialise_lattice(L);
    double beta = 100;
    int num_iterations = 10000;

    // Header of file
    ofs << "#t";
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            ofs << "#" << "i" << i << "j" << j;
        }
    }
    ofs << "\n";

    std::cout << "Starting simulation" << std::endl;
    for (int i = 0; i < num_iterations - 1000; ++i) {
        metropolis_algorithm(lattice, beta, J, H);
        writer(i, lattice, ofs);
        // std::cout << "step " << i << "\n";
    }
    for (int i = num_iterations - 1000; i < num_iterations; ++i) {
        writer(i, lattice, ofs);
        metropolis_algorithm(lattice, beta, J, H);
    }
    std::cout << "Simulation ended" << std::endl;
    ofs.close();
};
// TODO: add dimension parameter
// K(i,j) = i * L + j, always try to do row iteration
ivec initialise_lattice(int L) {
    ivec lattice;
    // generate random int array of 0 and 1
    lattice = rng::random_int_array(rng::engine, L * L);

    // (0,1)->(-1,1)
    for (auto &i : lattice) {
        i = (2 * i) - 1;
    }

    return lattice;
};

double compute_energy_at_site(int i, int j, ivec lattice) {
    double energy;
    int L = sqrt(lattice.size());
    int sign_change = lattice[i * L + j];

    // This is done using screw periodic boundary conditions.
    // TODO: Implement PBC and other types of boundary conditions
    int neighbours_sum =
        lattice[((i + 1) % L) * L + j] + lattice[((i + L - 1) % L) * L + j] +
        lattice[i * L + ((j + 1) % L)] + lattice[i * L + ((j + L - 1) % L)];

    energy = sign_change * neighbours_sum;
    return energy;
};

double compute_total_energy(ivec lattice, double J, double H) {
    double site_energy = 0;
    double magnetic_energy = 0;
    int L = sqrt(lattice.size());
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            site_energy += compute_energy_at_site(i, j, lattice);
            magnetic_energy += lattice[i * L + j];
        };
    };
    // Divide by 2 to account for over counting
    double total_energy = -J * (site_energy / 2) - H * magnetic_energy;
    return total_energy;
};

void metropolis_algorithm(ivec &lattice, double beta, double J, double H) {
    int L = sqrt(lattice.size());
    ivec indices = rng::random_int_array(rng::engine, 2, 0, L - 1);

    int i = indices[0];
    int j = indices[1];

    // Flip spin of site
    ivec new_lattice = lattice;
    new_lattice[i * L + j] = -new_lattice[i * L + j];

    // Calculate energy difference, by brute force
    // Calculate energy in a smarter way?
    double energy_diff = compute_total_energy(new_lattice, J, H) -
                         compute_total_energy(lattice, J, H);

    double r = rng::random_real_number(rng::engine);

    //  If r< e^-\beta \Delta E flip spin, ie update lattice
    double probability = exp(-beta * energy_diff);

    if (r < probability) {
        lattice[i * L + j] = -lattice[i * L + j];
    }
}

void writer(int iteration, ivec lattice, std::ofstream &file) {
    file << iteration;
    for (auto &i : lattice) {
        file << "," << i;
    }
    file << "\n";
}
