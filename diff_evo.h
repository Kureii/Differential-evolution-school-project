/*
 * Created by kureii on 11/15/24.
 */

#pragma once

#include <random>
#include <vector>

#include "structures.h"

class DiffEvo {
 public:
  DiffEvo();

  void GenerateInitPopulation(size_t dimensions, size_t populationSize);

  void SetInitPopulation(std::vector<de_t> initial_population);
  [[nodiscard]] const std::vector<de_t>& GetInitPopulation() const;

  void AddMinLimits(const std::vector<double>& min_limits);
  void AddMinLimits(double min_limits);
  void AddMaxLimits(const std::vector<double>& max_limits);
  void AddMaxLimits(double max_limits);

  template <typename FUN>
    requires std::invocable<FUN, de_t*>
  std::vector<de_t> Rand1(size_t iterations, double mutation_rate,
                          double cross_rate, FUN fun);

  template <typename FUN>
    requires std::invocable<FUN, de_t*>
  std::vector<de_t> Best1(size_t iterations, double mutation_rate,
                          double cross_rate, FUN fun);

  template <typename FUN>
    requires std::invocable<FUN, de_t *>
  std::vector<de_t> jDE(size_t iterations, double mutation_rate,
                        double cross_rate, double tau1, double tau2, FUN fun);

 private:
  std::vector<de_t> initial_population_;
  std::random_device rd_;
  std::mt19937 gen_;
  std::uniform_real_distribution<> distrib_;

  bool OptimizeInitTest(std::string& err);
};

#include "diff_evo.tpp"