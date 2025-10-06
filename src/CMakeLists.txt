# executables
add_executable(ising ising.cpp)

# make the executable in build directory
set_target_properties(ising PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)


# add them to global
target_include_directories(ising PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
