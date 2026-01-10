#ifndef MODEL_H
#define MODEL_H
#include "lattice/lattice.h"
#include "models/modelBase.h"
#include <iostream>
#include <string>
#include <vector>

using ivec = std::vector<int>;

template <typename T>
class Model: public ModelBase{
    protected:
        double beta,J, temp;
        // file object


    public: 
        //TODO manage the bullshit with virutal functions and destructors,
        // copies, etc...
        Model(double beta, double J, double temp)
            : beta(beta), J(J), temp(temp){}

        virtual T compute_total_spin_at_site(int index);
        virtual T compute_spin_neighbours_term(int index);
        virtual T compute_spin_neighbours_term(ivec indices);
        virtual T compute_spin_magnetic_term(int index);
        virtual T compute_spin_magnetic_term(ivec indices);
        void change_spin_randomly(int index) override;
        void change_spin_randomly(ivec indices) override;

        double compute_magnetisation() override;
        virtual T compute_total_energy();
        virtual T compute_energy_diff();

        // Class specific methods
        double compute_average_energy(std::vector<T> energies);
        double compute_heat_capacity(std::vector<T> energies);
        double compute_single_partition_function() override;

        void write_magnetisation() override;
        void write_energy() override;

};

#endif
