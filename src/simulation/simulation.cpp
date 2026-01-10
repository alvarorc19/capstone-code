#include "simulation/simulation.h"
#include "simulation/simulationParameters.h"
#include "simulation/simulationParameters.h"
#include "models/ising.h"
#include "models/potts.h"
#include "models/xy.h"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <hdf5.h>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>
#include <toml++/toml.hpp>

void Simulation::parse_parameters(std::string project_folder, std::string model_type) {
    toml::table tbl;

    // careful with this, has not been tested
    try
    {
        tbl = toml::parse_file(project_folder + "config.toml");
        std::cout << tbl << "\n";
    }
    catch (const toml::parse_error& err)
    {
        std::cout << "Make sure the parameters are in a file named config.toml" << "\n";
        std::cerr << "Parsing failed:\n" << err << "\n";
        return 2;
    }

    // Assign all the parameters to the struct
    parameters.L = static_cast<int>(tbl["physical_settings"]["L"].value_or<int64_t>(0));
    parameters.dim = static_cast<int>(tbl["physical_settings"]["dimension"].value_or<int64_t>(0));
    parameters.T = tbl["physical_settings"]["temperature"].value_or<double>(0.0);
    parameters.J = tbl["physical_settings"]["J"].value_or<double>(0.0);
    parameters.potts_q = static_cast<int>(tbl["physical_settings"]["potts_q"].value_or<int64_t>(0));
    parameters.model_type = model_type;
    parameters.time_steps = static_cast<int>(tbl["simulation_settings"]["time_steps"].value_or<int64_t>(0));
    parameters.recording_steps = static_cast<int>(tbl["simulation_settings"]["recording_steps"].value_or<int64_t>(0));
    parameters.record_magnetisation = tbl["simulation_settings"]["record_magnetisation"].value_or<bool>(false);
    parameters.record_energy = tbl["simulation_settings"]["record_energy"].value_or<bool>(false);
    parameters.record_susceptibility = tbl["simulation_settings"]["record_susceptibility"].value_or<bool>(false);
    parameters.record_specific_heat = tbl["simulation_settings"]["record_specific_heat"].value_or<bool>(false);
    parameters.record_correlation_length = tbl["simulation_settings"]["record_correlation_length"].value_or<bool>(false);
    parameters.record_correlation_function = tbl["simulation_settings"]["record_correlation_function"].value_or<bool>(false);}

void Simulation::run() {

    initialise_model();

    for(int i = 0; i < parameters.time_steps - parameters.recording_steps; i++){
        do_step();
        time_step++;
    }

    initialise_writing();

    for(int i = 0; i < parameters.recording_steps; i++){
        do_step();
        write();
        time_step++;
    }
}

void Simulation::initialise_model() {
    double beta = 1 / parameters.T;
    if (parameters.model_type == "ising"){
        model = std::make_unique<IsingModel>(
            beta,
            parameters.J,
            parameters.H,
            parameters.dim,
            parameters.L
        );
    }
    else if (parameters.model_type == "potts"){
        model = std::make_unique<PottsModel>(
            beta,
            parameters.J,
            // Figure this thing out
            parameters.vec_H,
            parameters.H,
            parameters.dim,
            parameters.L,
            parameters.potts_q
        );
    }
    else if (parameters.model_type == "xy"){
        model = std::make_unique<XYModel>(
            beta,
            parameters.J,
            // Figure this thing out
            parameters.vec_H,
            parameters.H,
            parameters.dim,
            parameters.L,
        );
    }
    else {
        throw std::invalid_argument("Invalid model_type: '" + parameters.model_type + "'. Must be 'ising', 'potts', or 'xy'.");
    }
}

void Simulation::initialise_writing() {
    int N = parameters.L * parameters.dim;
    std::string filename = folder_path + "results.h5";
    this-> file = HighFive::File(filename);

    std::vector<size_t> current_dims = {0, N}; 
    std::vector<size_t> max_dims = {HighFive::DataSpace::UNLIMITED, N}; 
    
    // Define chunks, how data is stored
    std::vector<size_t> chunk_dims = {parameters.recording_steps, N};

    HighFive::DataSpace lattice_space(current_dims, max_dims);
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(chunk_dims));

    parameters.lattice_dataset = file.createDataSet<int>(
        "lattice",
        lattice_space,
        props
    );

}

