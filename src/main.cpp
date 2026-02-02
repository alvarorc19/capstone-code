#include "simulation/simulation.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <memory>
#include <toml++/toml.hpp>
#include <omp.h>
#include <cstdio>

/**
 * @brief Main function to run the simulation
 * 
 * @param argc Number of arguments
 * @param argv Array of arguments
 * @return int Exit code
 */
int main(int argc, char *argv[]){
    std::string project_folder;
    std::string running_model;
    // Has the output and error to txt files
    // (void)std::freopen("output.txt", "w", stdout);
    // (void)std::freopen("error.txt", "w", stderr);
    int num_threads = 4;

        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];

            if (arg == "-p" && i + 1 < argc) {
                project_folder = argv[i + 1];
                ++i; // Skip next argument since it's the folder name
            }
            else if (arg == "-m" && i + 1 < argc){
                running_model = argv[i + 1];
                ++i;
            }
            else if (arg == "-j" && i + 1 < argc){
                num_threads = std::stoi(argv[i + 1]);
                ++i;
            }
        }

        if (!project_folder.empty()) {
            // Assuming it's being executed from build folder
            std::filesystem::path project_path = std::filesystem::current_path().parent_path() / std::filesystem::path(project_folder);
            std::cout << "Using project folder: " << project_path << "\n";
        }     

    
    try{ 
        omp_set_num_threads(num_threads);
        std::filesystem::path project_path = std::filesystem::current_path().parent_path() / std::filesystem::path(project_folder);
        std::vector<std::filesystem::path> directories;

        // Separate this for easier parallelisation implementation
        // Iterate over the different directories
        int i = 1;
        for (auto const& dir_entry : std::filesystem::directory_iterator{project_path}){ 
            if (dir_entry.is_directory()) {
                std::cout <<"Path of combination "<< i << ": " << dir_entry.path().filename() << '\n';
                directories.push_back(dir_entry.path());
                i++;
            }
        }
        std::cout << "Number of parameters = " << directories.size() << std::endl;
        /* #pragma omp parallel for schedule(dynamic) */
        #pragma omp parallel for
        for (int i = 0; i < directories.size(); i++){
            // out commands need to be wrapped around this
            #pragma omp critical
            {
            std::cout << "Starting parameter combination number " << i +1 << " out of " <<
                directories.size() << std::endl;
            }
            std::unique_ptr<Simulation> sim = std::make_unique<Simulation>();
            // Initialise in critical since multithreading messes up parsing and writing
            #pragma omp critical
            {
                sim->parse_parameters(directories[i], running_model);
                sim->initialise_writing();
            }
            sim->run();

            #pragma omp critical 
            {
            std::cout << "\n" <<"Ended parameter combination number " << i + 1 << " out of " << directories.size()<< std::endl;
            }
        }
        std::cout << "\n" << "Ended the simulation!" << std::endl;
        return 0;
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Configuration error: " << e.what() << std::endl;
        return 4;
    }
}
