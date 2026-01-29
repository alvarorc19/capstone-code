#ifndef SIMULATION_H
#define SIMULATION_H
#include "simulation/simulationParameters.h"
#include "simulation/observables.h"
#include "models/modelBase.h"
#include <iostream>
#include <vector>
#include <filesystem>
#include <memory>
#include <string>
#include <cmath>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>
#include <toml++/toml.hpp>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

/**
 * @class Simulation
 * @brief The Simulation class handles the overall simulation process,
 * including parameter parsing, model initialization, running the simulation,
 * and writing outputs.
 */
class Simulation {
    private:
        SimulationParameters parameters;/**< a SimulationParameters object that stores all the parameters needed for the simulation */ 
        Observables observables;/**< an Observables object that stores measurements during the simulation */
        std::unique_ptr<HighFive::File> file;/**< a unique_ptr element that points to the HDF5 output file */
        int time_step = 0;
        std::unique_ptr<ModelBase> model;/**< a unique_ptr element that points to any of the models created with the Model class or its parent ModelBase */

        /**
         * @brief Extract an integer vector from a toml array
         * 
         * @param array array inherited from toml obtained from the file
         * @return ivec data of toml array in int vector
         */
        ivec extract_int_vector_toml(const toml::array& array);

        /**
         * @brief Extract a double vector from a toml array
         * 
         * @param array array inherited from toml obtained from the file
         * @return dvec data of toml array in double vector
         */
        dvec extract_double_vector_toml(const toml::array& array);

    public:
        /**
         * @brief Construct a new Simulation object
         * 
         */
        Simulation() =default;

        /**
         * @brief Parse the parameters from the toml file
         * 
         * @param project_folder_path std::filesystem::path to the project folder
         * @param model_type std::string specifying the model type
         */
        void parse_parameters(std::filesystem::path project_folder_path, std::string model_type);

        /**
         * @brief Selects the model to be used
         * 
         * @detail Given a model parsed from the terminal it assigns the ModelBase pointer
         * to the correpsondent model class to execute the simulation
         * 
         */
        void initialise_model();

        /**
         * @brief Initialises writing objects
         * 
         * @detail Given the file pointer it creates the h5 file corresponding to that run.
         * It also creates the lattice dataset and dataspace as well as reserving space and defining the chunk size
         * for later writing.
         * 
         */
        void initialise_writing();
        
        /**
         * @brief It executes the simulation, a full run corresponding to the time steps
         * 
         */
        void run();

        /**
         * @brief executes a single Metropolis step
         * 
         */
        void do_metropolis_step();

        /**
         * @brief Executes a full sweep of Metropolis steps while recording observables
         * 
         * @detail TODO: explain algorithm
         * 
         * @param time_step int current time step
         * @return int returns the updated time step after the sweep
         */
        int do_metropolis_recording_sweep(int time_step);

        //TODO
        // void do_cluster_step();

        /**
         * @brief  Writes the lattice in the dataset at a given time step
         * 
         * @detail It resizes the dataset to accomodate the new time step and writes the lattice into the dataset
         * stored in the SimulationParameters struct.
         * 
         * @param time int time step
         */
        void write_lattice(int time);

        /**
         * @brief writes the observables into the h5 file after the run from the Observables object
         * 
         */
        void write_observables();

        /**
         * @brief Updates the observables stored in the Observables object
         * 
         */
        void update_observables();
};

#endif
