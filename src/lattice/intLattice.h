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
        ivec generate_lattice(int L, int dim, int q);

    public:
        IntLattice(int L , int dim, int q); // Done
        IntLattice(int L , int dim); // Done

        int get_lattice_site(ivec indices);
        int get_lattice_site(int index);

        ivec get_neighbours_array(ivec indices);
        ivec get_neighbours_array(int index);

};

#endif
