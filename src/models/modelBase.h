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
        virtual void change_spin(int random_index, double spin) = 0;
        virtual std::vector<int> cluster_flip_neighbours(int index, double direction, double angle_flip) = 0; 
        virtual void write_lattice(HighFive::DataSet * lattice_set, int time) = 0;
};

#endif
