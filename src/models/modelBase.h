#ifndef MODELBASE_H
#define MODELBASE_H
#include <iostream>

class ModelBase {
    public:
        virtual ~ModelBase() = default;
        virtual double compute_single_partition_function() = 0;
        virtual void revert_state_of_system() = 0;
        virtual double compute_total_energy() = 0;
        virtual double compute_magnetisation() = 0;
        virtual double compute_energy_diff_flip() = 0;
        virtual double compute_spin_magnetic_term(int dim) = 0;
        virtual void change_spin(int random_index, double spin) = 0;
        virtual void flip_spin(int index, double angle) = 0;
        virtual void cluster_flip_neighbours(int index, double direction, ivec& cluster_stack, int& spins_flipped, std::vector<uint8_t> & visited, int lattice_dim) = 0;
        virtual void write_lattice(HighFive::DataSet * lattice_set, hsize_t time, hsize_t N) = 0;
};

#endif
