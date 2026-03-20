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
    dvec rg_x_spins;
    dvec rg_y_spins;

  public:
    // XYModel(double temp, double J, dvec vec_H, int H, int dim, int L)
        // : Model<double>(temp, J), H(H),vec_H(vec_H) {
    XYModel(double temp, double J, int dim, int L)
        : Model<double>(temp, J) {
        this->lattice_obj = std::make_unique<DoubleLattice>(L,dim);
    }
        double compute_spin_neighbours_term(int index) override;
        double compute_spin_neighbours_term(ivec indices) override;
        double compute_spin_magnetic_term(int dim) override;
        double compute_spin_magnetic_term() override {return 0;};
        double compute_total_energy() override;
        double compute_magnetisation() override;
        double compute_energy_diff_flip() override;
        void cluster_flip_neighbours(int index, double direction, ivec& cluster_stack, int& spins_flipped, std::vector<uint8_t> & visited, int lattice_dim) override;
        void change_spin_randomly(ivec indices) override;
        void change_spin(int index, double spin) override;
        void flip_spin(int index, double angle) override;
        void compute_reduced_lattice(int b, size_t N, size_t L) override;
        double compute_rg_spin_magnetic_term(int dim) override;
        double compute_rg_energy(int b, size_t N, ivec & neigh_table, int dim) override;
        ivec calculate_reduced_neighbours_table(int L, int dim, int b) override;
};

#endif
