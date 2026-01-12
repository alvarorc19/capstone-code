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
       dvec generate_lattice(int L, int dim) override;

    public:
        DoubleLattice(const int L, const int dim);

        double& get_lattice_site(ivec indices) override;
        double& get_lattice_site(int index) override;

        dvec get_neighbours_array(ivec indices) override;
        dvec get_neighbours_array(int index) override;

};

#endif
