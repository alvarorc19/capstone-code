#include "lattice/lattice.h"
#include "random.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

template<class T>
int Lattice<T>::get_1d_index(ivec indices) {
    assert((lattice_dim == indices.size()) &&
           "Indices provided need to be same size as lattice");

    int index = 0;

    for (int i = 0; i < lattice_dim; i++) {
        assert((indices[i] < lattice_length) && "Indices must be in the lattice");
        index += indices[i] * std::pow(lattice_length, lattice_dim - i - 1);
    }

    return index;
}


template<class T>
ivec Lattice<T>::calculate_neighbours_table(int L, int dim) {
    ivec neighbours_array(std::pow(L, dim) * 2 * dim);
    if (dim == 2) {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {

                int index = i * L + j;

                // i
                int i_up = (i + 1 == L) ? 0 : i + 1;
                int i_down = (i - 1 == -1) ? L - 1 : i - 1;
                // j
                int j_up = (j + 1 == L) ? 0 : j + 1;
                int j_down = (j - 1 == -1) ? L - 1 : j - 1;

                // Up
                neighbours_array[4 * index + 0] = i_up * L + j;
                // Down
                neighbours_array[4 * index + 1] = i_down * L + j;
                // Right
                neighbours_array[4 * index + 2] = j_up * L + i;
                // Left
                neighbours_array[4 * index + 3] = j_down * L + i;
            }
        }

    } else if (dim == 3) {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                for (int k = 0; k < L; k++) {

                    int index = i * L * L + j * L + k;

                    // i
                    int i_up = (i + 1 == L) ? 0 : i + 1;
                    int i_down = (i - 1 == -1) ? L - 1 : i - 1;
                    // j
                    int j_up = (j + 1 == L) ? 0 : j + 1;
                    int j_down = (j - 1 == -1) ? L - 1 : j - 1;
                    // k
                    int k_up = (k + 1 == L) ? 0 : k + 1;
                    int k_down = (k - 1 == -1) ? L - 1 : k - 1;

                    // i
                    neighbours_array[6 * index + 0] = i_up * L * L + j * L + k;
                    neighbours_array[6 * index + 1] = i_down * L * L + j * L + k;
                    // j
                    neighbours_array[6 * index + 2] = i * L * L + j_up * L + k;
                    neighbours_array[6 * index + 3] = i * L * L + j_down * L + k;
                    // k
                    neighbours_array[6 * index + 4] = i * L * L + j * L + k_up;
                    neighbours_array[6 * index + 5] = i * L * L + j * L + k_down;
                }
            }
        }

    } else {
        // TODO Look at obsidian to know how to further implement it
        std::cout << "Further dimensions not implemented yet" << "\n";
        exit(2);
    }
    return neighbours_array;
}

