#ifndef MODEL_H
#define MODEL_H
#include "lattice/lattice.h"
#include "models/modelBase.h"
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>


using ivec = std::vector<int>;

template <typename T>
class Model: public ModelBase{
    protected:
        double beta,J, temp;
        std::unique_ptr<Lattice<T>> lattice_obj;
        std::vector<T> old_lattice;


    public: 
        //TODO manage the bullshit with virutal functions and destructors,
        // copies, etc...
        virtual ~Model() =default;
        Model(double temp, double J)
            : ModelBase(),beta(1/temp), J(J), temp(temp){}

        // Class specific methods
        double compute_average_energy(std::vector<T> energies);
        double compute_heat_capacity(std::vector<T> energies);
        double compute_single_partition_function() override;
        void revert_state_of_system() override {lattice_obj->get_lattice()=old_lattice;}
        void write_lattice(HighFive::DataSet& lattice_set, int time) override;

        // Virtual functions
        virtual T compute_spin_neighbours_term(int index) = 0;
        virtual T compute_spin_neighbours_term(ivec indices) = 0;
        virtual T compute_spin_magnetic_term(int dim) = 0;
        virtual T compute_spin_magnetic_term() = 0;
        virtual void change_spin_randomly(ivec indices) = 0;

};

#include "models/model.tpp"

#endif
