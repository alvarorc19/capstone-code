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

  public:
    XYModel(double temp, double J, dvec vec_H, int H, int dim, int L)
        : Model<double>(temp, J), H(H),vec_H(vec_H) {
        this->lattice_obj = std::make_unique<DoubleLattice>(L,dim);
    }
        double compute_spin_neighbours_term(int index) override;
        double compute_spin_neighbours_term(ivec indices) override;
        double compute_spin_magnetic_term(int dim) override;
        double compute_spin_magnetic_term() override {return 0;};
        double compute_total_energy() override;
        double compute_magnetisation() override;
        void cluster_flip_neighbours(int index, double direction, int& new_spins_flipped, double angle_flip) override;
        void change_spin_randomly(ivec indices) override;
    // TODO make it general
        void change_spin(int index, double spin) override;
        double compute_energy_diff_flip() override;
};

#endif
