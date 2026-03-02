#ifndef POTTSMODEL_H
#define POTTSMODEL_H
#include "lattice/intLattice.h"
#include "models/model.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using ivec = std::vector<int>;
using dvec = std::vector<double>;

class PottsModel : public Model<int> {
  private:
    dvec vec_H;
    double H;
    int q;

    int delta_function(int a, int b);

  public:
    // PottsModel(double temp, double J, dvec vec_H, int H, int dim, int L, int q)
    //     : Model<int>(temp, J), H(H),vec_H(vec_H), q(q) {
    PottsModel(double temp, double J, int dim, int L, int q)
        : Model<int>(temp, J), q(q) {
        this->lattice_obj = std::make_unique<IntLattice>(L,dim,q);
    }

        int compute_spin_neighbours_term(int index) override;
        int compute_spin_neighbours_term(ivec indices) override;
        // int compute_spin_magnetic_term(int dim) override;
        double compute_spin_magnetic_term(int dim) override{return 0.;};
        int compute_spin_magnetic_term() override{return 0;};
        double compute_total_energy() override;
        double compute_magnetisation() override;
        void change_spin_randomly(ivec indices) override;
        double compute_energy_diff_flip() override;
};

#endif
