/*
 * Created by kureii on 11/15/24.
 */

#include "diff_evo.h"

#include <algorithm>
#include <cfloat>
#include <stdexcept>

DiffEvo::DiffEvo() : gen_(rd_()), distrib_(0.0, 1.0) {}

const std::vector<de_t> &DiffEvo::GetInitPopulation() const {
  return initial_population_;
}

void DiffEvo::GenerateInitPopulation(size_t dimensions, size_t populationSize) {
  std::vector<de_t> population(populationSize);
  std::ranges::for_each(population, [&dimensions, this](auto &v) {
    de_t des;
    des.cost = DBL_MAX;
    des.weight.resize(dimensions);
    for (size_t i = 0; i < dimensions; ++i) {
      des.weight[i] = distrib_(gen_);
    }
    v = des;
  });
  initial_population_ = population;
}

void DiffEvo::SetInitPopulation(std::vector<de_t> initial_population) {
  initial_population_ = std::move(initial_population);
}

void DiffEvo::AddMinLimits(const std::vector<double> &min_limits) {
  if (min_limits.empty()) {
    throw std::invalid_argument("min_limits cannot be empty");
  }
  if (initial_population_.empty()) {
    throw std::invalid_argument(
        R"(initial_population must be initialized, use "DiffEvo::GenerateInitPopulation" or "DiffEvo::SetInitPopulation")");
  }
  if (min_limits.size() != initial_population_[0].weight.size()) {
    throw std::invalid_argument(
        "initial_population size does not match min_limits");
  }

  std::ranges::for_each(initial_population_,
                        [&min_limits](auto &v) { v.min_limits = min_limits; });
}

void DiffEvo::AddMinLimits(double min_limits) {
  if (initial_population_.empty()) {
    throw std::invalid_argument(
        R"(initial_population must be initialized, use "DiffEvo::GenerateInitPopulation" or "DiffEvo::SetInitPopulation")");
  }

  std::ranges::for_each(initial_population_, [&min_limits](auto &individual) {
    individual.min_limits.resize(individual.weight.size());
    std::ranges::for_each(individual.min_limits,
                          [&min_limits](double &value) { value = min_limits; });
  });
}

void DiffEvo::AddMaxLimits(const std::vector<double> &max_limits) {
  if (max_limits.empty()) {
    throw std::invalid_argument("max_limits cannot be empty");
  }
  if (initial_population_.empty()) {
    throw std::invalid_argument(
        R"(initial_population must be initialized, use "DiffEvo::GenerateInitPopulation" or "DiffEvo::SetInitPopulation")");
  }
  if (max_limits.size() != initial_population_[0].weight.size()) {
    throw std::invalid_argument(
        "initial_population weight size does not match max_limits");
  }

  std::ranges::for_each(initial_population_,
                        [&max_limits](auto &v) { v.max_limits = max_limits; });
}

void DiffEvo::AddMaxLimits(double max_limits) {
  if (initial_population_.empty()) {
    throw std::invalid_argument(
        R"(initial_population must be initialized, use "DiffEvo::GenerateInitPopulation" or "DiffEvo::SetInitPopulation")");
  }

  std::ranges::for_each(initial_population_, [&max_limits](auto &individual) {
    individual.max_limits.resize(individual.weight.size());

    std::ranges::for_each(individual.max_limits,
                          [&max_limits](double &value) { value = max_limits; });
  });
}

bool DiffEvo::OptimizeInitTest(std::string &err) {
  if (initial_population_.empty()) {
    err =
        R"(initial_population must be initialized, use "DiffEvo::GenerateInitPopulation" or "DiffEvo::SetInitPopulation")";
    return false;
  }

  if (initial_population_[0].min_limits.empty()) {
    err = R"(min_limits are empty, use "DiffEvo::AddMinLimits")";
    return false;
  }

  if (initial_population_[0].max_limits.empty()) {
    err = R"(max_limits are empty, use "DiffEvo::AddMaxLimits")";
    return false;
  }

  for (size_t i = 0; i < initial_population_[0].min_limits.size(); ++i) {
    if (initial_population_[0].min_limits[i] >
        initial_population_[0].max_limits[i]) {
      err = std::format("Limit on index {} is incorrect, {} !< {}", i,
                        initial_population_[0].min_limits[i],
                        initial_population_[0].max_limits[i]);
      return false;
    }
  }

  return true;
}
