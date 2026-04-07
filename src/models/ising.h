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
    // double H;

  public:
    // IsingModel(double temp, double J, double H, const int dim, const int L)
    //     : Model<int>(temp, J), H(H) {
    IsingModel(double temp, double J, const int dim, const int L)
        : Model<int>(temp, J) {
        this->lattice_obj = std::make_unique<IntLattice>(L,dim);
    }

    int compute_spin_neighbours_term(int index) override;
    int compute_spin_neighbours_term(ivec indices) override;
    int compute_spin_magnetic_term() override;
    double compute_spin_magnetic_term(int dim) override {return 0.;}
    double compute_total_energy() override;
    double compute_magnetisation() override;
    void change_spin_randomly(ivec indices) override;
    double compute_energy_diff_flip() override;
};


#endif
