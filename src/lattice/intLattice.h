#ifndef INT_LATTICE_H
#define INT_LATTICE_H
#include "lattice/lattice.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using ivec = std::vector<int>;

// Do not forget to make the code more robust and think were to use pointers
class IntLattice : public Lattice<int> {

    private:
        int potts_q;

        // Specifically for Potts' model
        ivec generate_lattice(int L, int dim, int q); 
        ivec generate_lattice(int L, int dim) override;

    public:
        IntLattice(int L , int dim, int q);
        IntLattice(int L , int dim);

        int& get_lattice_site(ivec indices) override;
        int& get_lattice_site(int index) override;

        ivec get_neighbours_array(ivec indices) override;
        ivec get_neighbours_array(int index) override;

};

#endif
