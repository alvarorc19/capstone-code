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

/**
 * @class Model
 * @brief Base clase for different models
 * 
 * This class provides the core functionality for simulating statistical mechanical
 * models (e.g., Ising, Potts, XY models) on lattices. It handles energy computations,
 * partition functions, and state management.
 * 
 * 
 * @tparam T data type used for lattice sites
 */
template <typename T>
class Model: public ModelBase{
    protected:
        double beta,J, temp;
        std::unique_ptr<Lattice<T>> lattice_obj;
        std::vector<T> old_lattice;


    public: 
        //TODO manage the virutal functions and destructors,
        // copies, etc...

        /**
         * @brief Destroy the Model object
         */
        virtual ~Model() =default;

        /**
         * @brief Default Constructor
         * 
         * Assigns the temperature and coupling constant to their private variables
         * 
         * @param temp 
         * @param J 
         */
        Model(double temp, double J)
            : ModelBase(),beta(1/temp), J(J), temp(temp){}

        // Class specific methods
        /**
         * @brief Computes the average energy of the system
         * 
         * @param energies array of energies to be averaged
         * @return double value of average energy 
         */
        double compute_average_energy(std::vector<T> energies);

        /**
         * @brief Computes the heat capacity of the system
         * @details TODO include formula
         * 
         * @param energies array of energies
         * @return double 
         */
        double compute_heat_capacity(std::vector<T> energies);

        /**
         * @brief Computes the partition function for a single configuration of the system
         * @details TODO INCLUDE FORMULA
         * 
         * @return double 
         */

        double compute_single_partition_function() override;
        
        /**
         * @brief Given an old configuration it reverts the system to that state
         * 
         */

        void revert_state_of_system() override {lattice_obj->get_lattice()=old_lattice;}
        
        /**
         * @brief Given a DataSet it writes the current time step lattice configuration to it
         * @details After opening hte h5 file and the dataset, this function is called to write the current
         * lattice configuration to it at the correct time step index.
         * 
         * @param[out] lattice_set a HighFive dataset object where the object is to be built
         * @param[in] time 
         */

        void write_lattice(HighFive::DataSet* lattice_set, hsize_t time, hsize_t N) override;

        // Virtual functions
        /**
         * @brief Computes the spin term due to neighbouring spins for a given index
         * 
         * @return T total spiin
         */
        virtual T compute_spin_neighbours_term(int index) = 0;

        /**
         * @brief Computes the spin term due to neighbouring spins for a given set of indices
         * 
         * @return T total spiin
         */
        virtual T compute_spin_neighbours_term(ivec indices) = 0;

        /**
         * @brief computes the spin in the magnetic term of the hamiltonian for a given dimension
         * 
         * @param dim dimension along which to compute the magnetic term
         * @return double total spin
         */
        virtual double compute_spin_magnetic_term(int dim) = 0;

        /**
         * @brief computes the spin in the magnetic term of the hamiltonian for a single dimension model
         * 
         * @return T total spin
         */
        virtual T compute_spin_magnetic_term() = 0;

        /**
         * @brief Given a site lattice index, it randomly changes the spin at that site
         * 
         * @param indices array of indices where the spins are to be changed
         */
        virtual void change_spin_randomly(ivec indices) = 0;
    // TODO implement this
        void change_spin(int index, double spin) override {}
        void flip_spin(int index, double angle) override {}
        void cluster_flip_neighbours(int index, double direction, ivec& cluster_stack, int& spins_flipped, std::vector<uint8_t> & visited, int lattice_dim) override{}
        void compute_reduced_lattice(int b, size_t N, size_t L) override{}
        void compute_reduced_lattice_3d(int b, size_t N, size_t L) override{}
        double compute_rg_spin_magnetic_term(int dim) override{return 1.;}
        double compute_rg_energy(int b, size_t N, ivec & neigh_table, int dim) override{return 0.;}
        ivec calculate_reduced_neighbours_table(int L, int dim, int b) override{ivec a{0};return a;}

};

#include "models/model.tpp"

#endif
