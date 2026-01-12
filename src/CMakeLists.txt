# executables
add_library(capstone_code STATIC 
    simulation/simulation.cpp
    random.cpp
    
    lattice/doubleLattice.cpp
    lattice/intLattice.cpp
    
    models/ising.cpp
    models/potts.cpp
    models/xy.cpp
)

# add them to global
target_include_directories(capstone_code PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(capstone_code
    PUBLIC
    HighFive
)


add_executable(run_simulation main.cpp)

target_link_libraries(run_simulation PRIVATE capstone_code)

# make the executable in build directory
set_target_properties(run_simulation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
install(TARGETS run_simulation DESTINATION .)

