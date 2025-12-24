#ifndef LATTICE_H
#define LATTICE_H
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <hdf5.h>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

// Implement & and stuff asa well as const and maybe virtual functions
template <typename T>
class Lattice {
    protected:
        // Variables or objects
        int lattice_dim;
        int lattice_length;
        ivec neighbours_table;
        std::vector<T> lattice;
        HighFive::File output_file;
        HighFive::DataSet lattice_set;

        // // Methods
        // Generate the lattices and tables
        virtual std::vector<T> generate_lattice(int L, int dim);
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
        Lattice(int L, int dim): lattice_dim(dim), lattice_length(L) {} // Done
        // Getters
        int get_1d_index(ivec indices);
        int get_lattice_total_length() {return std::pow(lattice_length, lattice_dim);};
        int get_lattice_length() { return lattice_length;}
        int get_lattice_dim() { return lattice_dim;}
        virtual const std::vector<T>& get_lattice() const {return lattice;}
        virtual std::vector<T>& get_lattice(){return lattice;}

        //Type specific
        virtual T get_lattice_site(ivec indices);
        virtual T get_lattice_site(int index);

        virtual std::vector<T> get_neighbours_array(ivec indices);
        virtual std::vector<T> get_neighbours_array(int index);

        void writer();
};
#endif
