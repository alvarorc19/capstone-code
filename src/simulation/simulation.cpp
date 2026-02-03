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
#include <chrono>
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
    parameters.dim = static_cast<int>(tbl["physical_settings"]["dimension"].value_or<int64_t>(0));
    parameters.L = static_cast<int>(tbl["physical_settings"]["L"].value_or<int64_t>(0));
    parameters.N = std::pow(parameters.L, parameters.dim);
    parameters.T = tbl["physical_settings"]["temperature"].value_or<double>(0.0);
    parameters.beta = 1 / parameters.T;
    parameters.J = tbl["physical_settings"]["J"].value_or<double>(0.0);
    parameters.potts_q = static_cast<int>(tbl["physical_settings"]["potts_q"].value_or<int64_t>(0));
    parameters.H = tbl["physical_settings"]["H"].value_or<double>(0.0);

    parameters.vec_H = extract_double_vector_toml(*tbl["physical_settings"]["vec_H"].as_array());

    parameters.model_type = model_type;
    parameters.total_sweeps = static_cast<size_t>(tbl["simulation_settings"]["total_sweeps"].value_or<int64_t>(0));
    parameters.recording_sweeps = static_cast<size_t>(tbl["simulation_settings"]["recording_sweeps"].value_or<int64_t>(0));
    parameters.record_lattice = tbl["simulation_settings"]["record_lattice"].value_or<bool>(false);
    parameters.record_correlation_length = tbl["simulation_settings"]["record_correlation_length"].value_or<bool>(false);
    parameters.record_correlation_function = tbl["simulation_settings"]["record_correlation_function"].value_or<bool>(false);}

void Simulation::run() {

    visited.assign(parameters.N, 0);
    cluster_stack.reserve(parameters.N);
    // Start timer
    const auto start{std::chrono::steady_clock::now()};

    do_cluster_sweep();

    // End timer
    const auto finish{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{finish - start};
    time_step++;

    #pragma omp critical
    {
        std::cout << "One sweep took " << elapsed_seconds.count() << " s\n";
        std::cout << "Expected execution time (not accounting for writing) =  " << 
            elapsed_seconds.count() * parameters.total_sweeps << " s " << std::endl;
    }

    for(int i = 1; i < parameters.total_sweeps - parameters.recording_sweeps; i++){
        // do_metropolis_sweep();
        do_cluster_sweep();
        time_step++;
        // #pragma omp critical
        // {
        //     std::cout << "Thermalisation step temp = " << parameters.T <<" = " << i <<"\n";
        // }
    }
    #pragma omp critical
    {
        std::cout << "Finished thermalisation, starting writing" << "\n";
    }

    for(int i = 0; i < parameters.recording_sweeps; i++){
        //time_step = do_metropolis_recording_sweep(time_step);
        //write_lattice(i);
        #pragma omp critical
        {
            std::cout << "step " << i << " out of " << parameters.recording_sweeps << "\n";
        }
        if (parameters.record_lattice and (parameters.T == 0.3 or parameters.T == 0.868421052631579 or parameters.T == 0.9315789473684211 or parameters.T == 1.5)) {
            for(int i = 0; i < parameters.recording_sweeps; i++){
                // do_metropolis_recording_sweep();
                // time_step++;
                // write_lattice(i);

                // Recording sweep records every step instead of every sweep (MCS)
                // do_cluster_recording_sweep();
                do_cluster_sweep();
                update_observables();
                write_lattice(i);
                time_step++;
                // #pragma omp critical
                // {
                // std::cout << "step " << i << " out of " << parameters.recording_steps << "\n";
                // }
            }
        }
        else if (not parameters.record_lattice) {
            for(int i = 0; i < parameters.recording_sweeps; i++){
                // do_metropolis_recording_sweep();
                
                // Recording sweep records every step instead of every sweep (MCS)
                // do_cluster_recording_sweep();
                do_cluster_sweep();
                update_observables();
                time_step++;
                // #pragma omp critical
                // {
                // std::cout << "step " << i << " out of " << parameters.recording_steps << "\n";
                // }
            }
        }
        else {
            for(int i = 0; i < parameters.recording_sweeps; i++){
                // do_metropolis_recording_sweep();
                
                // Recording sweep records every step instead of every sweep (MCS)
                // do_cluster_recording_sweep();
                do_cluster_sweep();
                update_observables();
                time_step++;

                if (spins_flipped > 1000 * parameters.N){
                    break;
                }
                // #pragma omp critical
                // {
                // std::cout << "step " << i << " out of " << parameters.recording_steps << "\n";
                // }
        }

        }
    }
    write_observables();
    file->flush();
    #pragma omp critical
    {
        std::cout << "Finished run for " << parameters.project_folder_path.filename() << "\n";
    }
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
    std::filesystem::path filename = parameters.project_folder_path / std::filesystem::path("results.h5");
    this->file = std::make_unique<HighFive::File>(filename.string(), HighFive::File::Truncate);

    if (parameters.record_lattice){

    std::vector<size_t> current_dims = {0, parameters.N}; 
    std::vector<size_t> max_dims = {HighFive::DataSpace::UNLIMITED, parameters.N}; 
    
    // Define chunks, how data is stored
    size_t chunk_rows = std::min(parameters.recording_sweeps, static_cast<size_t>(375));
    std::vector<hsize_t> chunk_dims = {chunk_rows, parameters.N};

    HighFive::DataSpace lattice_space(current_dims, max_dims);
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(chunk_dims));

    if (parameters.model_type == "xy") {
        parameters.lattice_set = std::make_unique<HighFive::DataSet>(
            file->createDataSet<double>(
                "lattice",
                lattice_space,
                props
            )
        );
    }
    else if (parameters.model_type == "ising" or parameters.model_type == "potts") {
        parameters.lattice_set = std::make_unique<HighFive::DataSet>(
            file->createDataSet<int>(
                "lattice",
                lattice_space,
                props
            )
        );
    }
    else { 
        throw std::invalid_argument("There was an oopsie when there shouldn't, I do not know how you got here");
    }
    }
    
    // Beware with the space reserved
    observables.energy_array.reserve(parameters.recording_sweeps * parameters.N);
    observables.magnetisation_array.reserve(parameters.recording_sweeps * parameters.N);
}

void Simulation::do_metropolis_step() {
    double energy_diff = model->compute_energy_diff_flip();
    double r = rng::random_real_number(rng::engine);
    double probability = exp(-parameters.beta * energy_diff);
    if (r > probability) {
        model->revert_state_of_system();
    }
}

/** 
* Here time is measured as MCS/site, Monte Carlo Steps per site.
* So 1 step is N flips for N particles
*/
void Simulation::do_metropolis_sweep(){
    for (int i = 0; i < parameters.N; i++){
        do_metropolis_step();
    }
}

void Simulation::do_metropolis_recording_sweep() {
    for (int i = 0; i < parameters.N; i++){
        do_metropolis_step();
        update_observables();
    }
}

/**
* Here we need to normalise the MCS steps things (or not)
* I can do this like this at the moment and leave it.
*
* TODO Normalise the MCS,
* 
* Total MCS = {\sum (size of each cluster)} / number particles
*/
void Simulation::do_cluster_step(){
    // Get the random stuff
    int random_index = rng::random_int_number(rng::engine, 0, parameters.N-1);
    double angle_r = rng::random_angle(rng::engine);
    double angle_flip = rng::random_angle(rng::engine);
    model->change_spin(random_index, angle_flip);
    
    // Make stack to store indices used
    cluster_stack.clear();
    cluster_stack.emplace_back(random_index);
    std::fill(visited.begin(), visited.end(), 0);
    cluster_stack.reserve(parameters.N);

    int current_index;
    spins_flipped++;

    while(!cluster_stack.empty()){
        current_index = cluster_stack.back();
        cluster_stack.pop_back();
        model->cluster_flip_neighbours(current_index, angle_r,angle_flip, cluster_stack, spins_flipped, visited);
        
    } 
}


void Simulation::do_cluster_sweep(){
    spins_flipped = 0;
        while (spins_flipped < parameters.N) {
            do_cluster_step();
        }
}

void Simulation::do_cluster_recording_sweep(){
    for (int i = 0; i < parameters.N; i++){
        do_cluster_step();
        update_observables();
        write_lattice(i+time_step);
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
    model->write_lattice(parameters.lattice_set.get(),time);
}

void Simulation::write_observables() {
    HighFive::DataSet magnetisation_dataset = file->createDataSet("observables/magnetisation", observables.magnetisation_array);
    HighFive::DataSet energy_dataset = file->createDataSet("observables/energy", observables.energy_array);
}
