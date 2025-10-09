#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

std::vector<int> initialise_lattice(int L);
double compute_energy_at_site(int x, int y, std::vector<int> lattice);
double compute_total_energy(std::vector<int> lattice);
void metropolis_algorithm(std::vector<int> &lattice, double beta);
void writer(int iteration, std::vector<int> lattice, std::ofstream &file);

int main(int argc, char *argv[]) {
    std::ofstream ofs("results.csv");

    int L = 100;
    std::vector<int> lattice;
    lattice = initialise_lattice(L);
    double beta = 10000;
    int num_iterations = 100000;

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
        metropolis_algorithm(lattice, beta);
    }
    for (int i = num_iterations - 1000; i < num_iterations; ++i) {
        writer(i, lattice, ofs);
        metropolis_algorithm(lattice, beta);
    }
    std::cout << "Simulation ended" << std::endl;
    ofs.close();
};
// TODO: add dimension parameter
// K(i,j) = i * L + j, always try to do row iteration
std::vector<int> initialise_lattice(int L) {
    std::vector<int> lattice;
    // generate random int array of 0 and 1
    lattice = rng::random_int_array(rng::engine, L * L);

    // (0,1)->(-1,1)
    for (auto &i : lattice) {
        i = (2 * i) - 1;
    }

    return lattice;
};

double compute_energy_at_site(int i, int j, std::vector<int> lattice) {
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

double compute_total_energy(std::vector<int> lattice) {
    double total_energy = 0;
    int L = sqrt(lattice.size());
    for (int i; i < L; ++i) {
        for (int j; j < L; ++j) {
            total_energy += compute_energy_at_site(i, j, lattice);
        };
    };
    // Divide by 2 to account for over counting
    return total_energy / 2;
};

void metropolis_algorithm(std::vector<int> &lattice, double beta) {
    int L = sqrt(lattice.size());
    std::vector<int> indices = rng::random_int_array(rng::engine, 2, 0, L - 1);

    int i = indices[0];
    int j = indices[1];

    // Flip spin of site
    std::vector<int> new_lattice = lattice;
    new_lattice[i * L + j] = -new_lattice[i * L + j];

    // Calculate energy difference, by brute force
    // Calculate energy in a smarter way?
    double energy_diff =
        compute_total_energy(new_lattice) - compute_total_energy(lattice);

    double r = rng::random_real_number(rng::engine);

    //  If r< e^-\beta \Delta E flip spin, ie update lattice
    double probability = exp(-beta * energy_diff);

    if (r < probability) {
        lattice[i * L + j] = -lattice[i * L + j];
    }
}

void writer(int iteration, std::vector<int> lattice, std::ofstream &file) {
    file << iteration;
    for (auto &i : lattice) {
        file << "," << i;
    }
    file << "\n";
}
