#include "simulation/simulation.h"
#include "simulation/simulationParameters.h"
#include "simulation/simulationParameters.h"
#include "models/ising.h"
#include "models/potts.h"
#include "models/xy.h"
#include "random.h"

#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <cmath>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>
#include <toml++/toml.hpp>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

void Simulation::parse_parameters(std::filesystem::path project_folder_path, std::string model_type) {
    toml::table tbl;
    // std::filesystem::path project_folder_path = std::filesystem::current_path().parent_path() / std::filesystem::path(project_folder);
    // careful with this, has not been tested
    try
    {
        std::filesystem::path config_path = project_folder_path / std::filesystem::path("config.toml");
        tbl = toml::parse_file(config_path.string());
        // std::cout << tbl << "\n";
    }
    catch (const toml::parse_error& err)
    {
        std::cout << "Make sure the parameters are in a file named config.toml" << "\n";
        std::cerr << "Parsing failed:\n" << err << "\n";
        // return 2;
    }

    // Assign the seed
    int seed = static_cast<int>(tbl["seed"].value_or<int64_t>(0));
    rng::update_seed(seed);

    // Assign all the parameters to the struct
    parameters.project_folder_path = project_folder_path;
    parameters.L = static_cast<int>(tbl["physical_settings"]["L"].value_or<int64_t>(0));
    parameters.dim = static_cast<int>(tbl["physical_settings"]["dimension"].value_or<int64_t>(0));
    parameters.T = tbl["physical_settings"]["temperature"].value_or<double>(0.0);
    parameters.beta = 1 / parameters.T;
    parameters.J = tbl["physical_settings"]["J"].value_or<double>(0.0);
    parameters.potts_q = static_cast<int>(tbl["physical_settings"]["potts_q"].value_or<int64_t>(0));
    parameters.H = tbl["physical_settings"]["H"].value_or<double>(0.0);

    parameters.vec_H = extract_double_vector_toml(*tbl["physical_settings"]["vec_H"].as_array());

    parameters.model_type = model_type;
    parameters.time_steps = static_cast<int>(tbl["simulation_settings"]["time_steps"].value_or<int64_t>(0));
    parameters.recording_steps = static_cast<size_t>(tbl["simulation_settings"]["recording_steps"].value_or<int64_t>(0));
    parameters.record_magnetisation = tbl["simulation_settings"]["record_magnetisation"].value_or<bool>(false);
    parameters.record_energy = tbl["simulation_settings"]["record_energy"].value_or<bool>(false);
    parameters.record_susceptibility = tbl["simulation_settings"]["record_susceptibility"].value_or<bool>(false);
    parameters.record_specific_heat = tbl["simulation_settings"]["record_specific_heat"].value_or<bool>(false);
    parameters.record_correlation_length = tbl["simulation_settings"]["record_correlation_length"].value_or<bool>(false);
    parameters.record_correlation_function = tbl["simulation_settings"]["record_correlation_function"].value_or<bool>(false);}

void Simulation::run() {

    initialise_model();

    for(int i = 0; i < parameters.time_steps - parameters.recording_steps; i++){
        do_metropolis_step();
        time_step++;
    }

    initialise_writing();

    for(int i = 0; i < parameters.recording_steps; i++){
        do_metropolis_step();
        write_lattice(time_step);
        update_observables();
        time_step++;
    }

    write_observables();
}

void Simulation::initialise_model() {
    if (parameters.model_type == "ising"){
        model = std::make_unique<IsingModel>(
            parameters.T,
            parameters.J,
            parameters.H,
            parameters.dim,
            parameters.L
        );
    }
    else if (parameters.model_type == "potts"){
        model = std::make_unique<PottsModel>(
            parameters.T,
            parameters.J,
            parameters.vec_H,
            parameters.H,
            parameters.dim,
            parameters.L,
            parameters.potts_q
        );
    }
    else if (parameters.model_type == "xy"){
        model = std::make_unique<XYModel>(
            parameters.T,
            parameters.J,
            parameters.vec_H,
            parameters.H,
            parameters.dim,
            parameters.L
        );
    }
    else {
        throw std::invalid_argument("Invalid model_type: '" + parameters.model_type + "'. Must be 'ising', 'potts', or 'xy'.");
    }
}

void Simulation::initialise_writing() {
    size_t N = std::pow(parameters.L, parameters.dim);
    std::filesystem::path filename = parameters.project_folder_path / std::filesystem::path("results.h5");
    this->file = std::make_unique<HighFive::File>(filename.string(), HighFive::File::Truncate);

    std::vector<size_t> current_dims = {0, N}; 
    std::vector<size_t> max_dims = {HighFive::DataSpace::UNLIMITED, N}; 
    
    // Define chunks, how data is stored
    size_t chunk_rows = std::min(parameters.recording_steps, static_cast<size_t>(375));
    std::vector<hsize_t> chunk_dims = {chunk_rows, N};

    HighFive::DataSpace lattice_space(current_dims, max_dims);
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(chunk_dims));

    if (parameters.model_type == "xy") {
        parameters.lattice_set = file->createDataSet<double>(
            "lattice",
            lattice_space,
            props
        );
    }
    else if (parameters.model_type == "ising" or parameters.model_type == "potts") {
        parameters.lattice_set = file->createDataSet<int>(
            "lattice",
            lattice_space,
            props
        );
    }
    else { 
        std::cout << "IDK bro" << "\n";
    }
    
    observables.energy_array.reserve(parameters.recording_steps);
    observables.magnetisation_array.reserve(parameters.recording_steps);
}

void Simulation::do_metropolis_step() {
    double energy_diff = model->compute_energy_diff_flip();
    double r = rng::random_real_number(rng::engine);
    double probability = exp(-parameters.beta * energy_diff);
    if (r > probability) {
        model->revert_state_of_system();
    }
}

ivec Simulation::extract_int_vector_toml(const toml::array & array){
    ivec int_vector;
    int_vector.reserve(array.size());

    for (const auto &i : array) {
        int value = i.as_integer()->get();
        int_vector.emplace_back(value);
    }
    return int_vector;
}

dvec Simulation::extract_double_vector_toml(const toml::array& array) {
    dvec double_vector;
    double_vector.reserve(array.size());
    
    for (const auto &i: array) {
        double value = i.as_floating_point()->get();
        double_vector.emplace_back(value);
    }

    return double_vector;
}

void Simulation::update_observables() {
    observables.energy_array.emplace_back(model->compute_total_energy());
    observables.magnetisation_array.emplace_back(model->compute_magnetisation());
}

void Simulation::write_lattice(int time) {
    model->write_lattice(parameters.lattice_set,time);
}

void Simulation::write_observables() {
    HighFive::DataSet magnetisation_dataset = file->createDataSet("observables/magnetisation", observables.magnetisation_array);
    HighFive::DataSet energy_dataset = file->createDataSet("observables/energy", observables.energy_array);
}
