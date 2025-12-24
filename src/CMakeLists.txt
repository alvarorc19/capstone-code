# executables
add_library(capstone_code STATIC 
    simulation/simulation.cpp 
    random.cpp
    
    lattice/doubleLattice.cpp
    lattice/intLattice.cpp
    lattice/lattice.cpp
    
    models/ising.cpp
    models/model.cpp
    models/potts.cpp
    models/xy.cpp
)

add_executable(run_simulation main.cpp)
# make the executable in build directory
set_target_properties(run_simulation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
install(TARGETS run_simulation DESTINATION .)

# add them to global
target_include_directories(run_simulation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
