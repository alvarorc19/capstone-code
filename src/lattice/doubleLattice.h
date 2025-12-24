#ifndef DOUBLE_LATTICE_H
#define DOUBLE_LATTICE_H
#include "lattice/lattice.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using dvec = std::vector<double>;
using ivec = std::vector<int>;

class DoubleLattice : public Lattice<double> {

    private:
       dvec generate_lattice(int L, int dim);

    public:
        DoubleLattice(const int L, const int dim);

        double get_lattice_site(ivec indices);
        double get_lattice_site(int index);

        std::vector<double> get_neighbours_array(ivec indices);
        std::vector<double> get_neighbours_array(int index);

};

#endif
