#ifndef ISINGMODEL_H
#define ISINGMODEL_H
#include "lattice/intLattice.h"
#include "models/model.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using ivec = std::vector<int>;

class IsingModel : public Model<int> {
  private:
    double H;
    IntLattice lattice_obj;

  public:
    IsingModel(double beta, double J, double H, const int dim, const int L)
        : Model(beta, J), H(H) {
        lattice_obj = IntLattice(dim, L)
    }

    int compute_total_spin_at_site(int index);
    int compute_spin_neighbours_term(int index);
    int compute_spin_neighbours_term(ivec indices);
    int compute_spin_magnetic_term(int index);
    int compute_spin_magnetic_term(ivec indices);
    void change_spin_randomly(int index);
    void change_spin_randomly(ivec indices);

    double compute_magnetisation();
    int compute_total_energy();
    int compute_energy_diff();
};


#endif
