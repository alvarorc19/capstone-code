#ifndef XYMODEL_H
#define XYMODEL_H
#include "lattice/doubleLattice.h"
#include "models/model.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

class XYModel : public Model<double> {
  private:
    dvec vec_H;
    double H;
    int q;
    DoubleLattice lattice_obj;

  public:
    XYModel(double beta, double J, dvec vec_H, int H, int dim, int L, int q)
        : Model(beta, J), H(H),vec_H(vec_H), q(q) {
        lattice_obj = DoubleLattice(dim, L)
    }
        double compute_total_spin_at_site(int index);
        double compute_spin_neighbours_term(int index);
        double compute_spin_neighbours_term(ivec indices);
        double compute_spin_magnetic_term(int index);
        double compute_spin_magnetic_term(ivec indices);
        void change_spin_randomly(int index);
        void change_spin_randomly(ivec indices);

        double compute_magnetisation();
        double compute_total_energy();
        double compute_energy_diff();
};

#endif
