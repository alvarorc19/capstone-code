#ifndef LATTICE_H
#define LATTICE_H
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

// Implement & and stuff asa well as const and maybe virtual functions
template <typename T>
class Lattice {
    protected:
        // Variables or objects
        const int lattice_dim;
        const int lattice_length;
        const int particle_num;
        ivec neighbours_table;
        std::vector<T> lattice;

        // // Methods
        // Generate the lattices and tables
        virtual std::vector<T> generate_lattice(int L, int dim) = 0;
        ivec calculate_neighbours_table(int L, int dim);

    public:
        // // Probably don't need this
        // // Rule of 5
        // Lattice() = delete;
        // Lattice(const Lattice &) = delete;
        // Lattice &operator=(const Lattice &other) = delete;
        // Lattice(Lattice &&) = delete;
        // Lattice &operator=(Lattice &&) = delete;
        virtual ~Lattice() =default;

        // Constructor
        Lattice(int L, int dim): lattice_dim(dim), lattice_length(L), particle_num(std::pow(L, dim)) {} // Done

        // Getters
        int get_1d_index(ivec indices);
        //TODO is this redundant??
        int get_lattice_total_length() {return std::pow(lattice_length, lattice_dim);};
        int get_lattice_length() { return lattice_length;}
        int get_lattice_dim() { return lattice_dim;}
        int get_particle_num() { return particle_num;}
        const std::vector<T>& get_lattice() const {return lattice;}
        std::vector<T>& get_lattice(){return lattice;}
        std::vector<int> get_neighbours_indices(ivec indices);
        std::vector<int> get_neighbours_indices(int index);

        //Type specific
        virtual T& get_lattice_site(ivec indices) = 0;
        virtual T& get_lattice_site(int index) = 0;

        virtual std::vector<T> get_neighbours_array(ivec indices) = 0;
        virtual std::vector<T> get_neighbours_array(int index) = 0;

};

#include "lattice/lattice.tpp"

#endif
