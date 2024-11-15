#pragma once

#include <algorithm>
#include <format>
#include <stdexcept>

template <typename FUN>
  requires std::invocable<FUN, de_t *>
std::vector<de_t> DiffEvo::Rand1( size_t iterations,
                                  double mutation_rate,
                                  double cross_rate, FUN fun) {

  if (std::string err_msg; !OptimizeInitTest(err_msg)) {
    throw std::runtime_error(err_msg);
  }

  std::vector<de_t> output;
  auto start_population = initial_population_;

  if (start_population.size() < 4) return output;

  const auto size = start_population.size();
  const auto dim = start_population[0].weight.size();

  std::ranges::for_each(start_population, [](auto &individual) {
    for (int i = 0; i < individual.weight.size(); ++i) {
        if (individual.weight[i] < individual.min_limits[i]) {
          individual.weight[i] = individual.min_limits[i];
        }
        if (individual.weight[i] > individual.max_limits[i]) {
          individual.weight[i] = individual.max_limits[i];
        }

    }
  });

  std::ranges::for_each(start_population, [&](auto &v) { fun(&v); });

  auto counter = start_population.size();
  std::vector<de_t> population(size);
  while (counter < iterations) {
    population.clear();
    while (population.size() < size) {
      de_t new_des;
      new_des.weight.resize(dim);
      new_des.min_limits = start_population[0].min_limits;
      new_des.max_limits = start_population[0].max_limits;

      std::vector<int> indices(start_population.size());
      std::ranges::iota(indices, 0);

      std::ranges::shuffle(indices, gen_);

      int x1_idx = indices[0];
      int x2_idx = indices[1];
      int x3_idx = indices[2];

      // mutate
      for (int i = 0; i < dim; ++i) {
        new_des.weight[i] =
            start_population[x1_idx].weight[i] +
            mutation_rate * (start_population[x2_idx].weight[i] -
                             start_population[x3_idx].weight[i]);
      }

      // crossing
      for (int i = 0; i < dim; ++i) {
        new_des.weight[i] = distrib_(gen_) > cross_rate
                                ? start_population[x1_idx].weight[i]
                                : new_des.weight[i];

        if (new_des.weight[i] < new_des.min_limits[i]) {
          new_des.weight[i] = new_des.min_limits[i];
        }
        if (new_des.weight[i] > new_des.max_limits[i]) {
          new_des.weight[i] = new_des.max_limits[i];
        }
      }

      fun(&new_des);
      population.push_back(new_des);
      counter++;
    }

    start_population = population;
    output.emplace_back(*std::ranges::min_element(
        start_population,
        [&](const auto &a, const auto &b) { return a.cost < b.cost; }));
  }

  return output;
}

template <typename FUN>
  requires std::invocable<FUN, de_t *>
std::vector<de_t> DiffEvo::Best1( size_t iterations,
                                  double mutation_rate,
                                  double cross_rate, FUN fun) {
  if (std::string err_msg; !OptimizeInitTest(err_msg)) {
    throw std::runtime_error(err_msg);
  }

  std::vector<de_t> output;
  auto start_population = initial_population_;

  if (start_population.size() < 4) return output;

  const auto size = start_population.size();
  const auto dim = start_population[0].weight.size();

  std::ranges::for_each(start_population, [](auto &individual) {
    for (int i = 0; i < individual.weight.size(); ++i) {
        if (individual.weight[i] < individual.min_limits[i]) {
          individual.weight[i] = individual.min_limits[i];
        }
        if (individual.weight[i] > individual.max_limits[i]) {
          individual.weight[i] = individual.max_limits[i];
        }

    }
  });

  std::ranges::for_each(start_population, [&](auto &v) { fun(&v); });

  auto counter = start_population.size();
  std::vector<de_t> population(size);
  while (counter < iterations) {
    population.clear();
    while (population.size() < size) {
      de_t new_des;
      new_des.weight.resize(dim);
      new_des.min_limits = start_population[0].min_limits;
      new_des.max_limits = start_population[0].max_limits;

      auto best = std::ranges::min_element(
          start_population,
          [&](const auto &a, const auto &b) { return a.cost < b.cost; });

      // gen indexes
      const auto best_idx = std::ranges::distance(start_population.begin(), best);
      std::vector<int> indices(start_population.size());
      std::iota(indices.begin(), indices.end(), 0);

      indices.erase(indices.begin() + best_idx);
      std::ranges::shuffle(indices, gen_);

      const int x2_idx = indices[0];
      const int x3_idx = indices[1];

      // mutation
      for (int i = 0; i < dim; ++i) {
        new_des.weight[i] =
            best->weight[i] +
            mutation_rate * (start_population[x2_idx].weight[i] -
                             start_population[x3_idx].weight[i]);
      }

      // crossing
      for (int i = 0; i < dim; ++i) {
        new_des.weight[i] =
            distrib_(gen_) > cross_rate ? best->weight[i] : new_des.weight[i];
        if (new_des.weight[i] < new_des.min_limits[i]) {
          new_des.weight[i] = new_des.min_limits[i];
        }
        if (new_des.weight[i] > new_des.max_limits[i]) {
          new_des.weight[i] = new_des.max_limits[i];
        }
      }
      fun(&new_des);
      population.push_back(new_des);
      counter++;
    }

    start_population = population;
    output.emplace_back(*std::ranges::min_element(
        start_population,
        [&](const auto &a, const auto &b) { return a.cost < b.cost; }));
  }

  return output;
}


template <typename FUN>
  requires std::invocable<FUN, de_t *>
std::vector<de_t> DiffEvo::jDE( size_t iterations,
                                  double mutation_rate,
                                  double cross_rate, double tau1, double tau2,
                                  FUN fun) {
  if (std::string err_msg; !OptimizeInitTest(err_msg)) {
    throw std::runtime_error(err_msg);
  }

  std::vector<de_t> output;
  auto start_population = initial_population_;
  std::uniform_real_distribution f_dist(0.1, 0.9);


  if (start_population.size() < 4) return output;

  const auto size = start_population.size();
  const auto dim = start_population[0].weight.size();

  std::ranges::for_each(start_population, [](auto &individual) {
    for (int i = 0; i < individual.weight.size(); ++i) {
        if (individual.weight[i] < individual.min_limits[i]) {
          individual.weight[i] = individual.min_limits[i];
        }
        if (individual.weight[i] > individual.max_limits[i]) {
          individual.weight[i] = individual.max_limits[i];
        }

    }
  });

  std::ranges::for_each(start_population, [&](auto &v) { fun(&v); });

  auto counter = start_population.size();
  std::vector<de_t> population(size);
  while (counter < iterations) {
    population.clear();
    while (population.size() < size) {
      de_t new_des;
      new_des.weight.resize(dim);
      new_des.min_limits = start_population[0].min_limits;
      new_des.max_limits = start_population[0].max_limits;

      if (distrib_(gen_) < tau1) {
        cross_rate = distrib_(gen_);
      }
      if (distrib_(gen_) < tau2) {
        mutation_rate = f_dist(gen_);
      }

      // gen indexes
      std::vector<int> indices(start_population.size());
      std::ranges::iota(indices, 0);

      std::ranges::shuffle(indices, gen_);

      int x1_idx = indices[0];
      int x2_idx = indices[1];
      int x3_idx = indices[2];

      // mutate
      for (int i = 0; i < dim; ++i) {
        new_des.weight[i] =
            start_population[x1_idx].weight[i] +
            mutation_rate * (start_population[x2_idx].weight[i] -
                             start_population[x3_idx].weight[i]);
      }

      // crossing
      for (int i = 0; i < dim; ++i) {
        new_des.weight[i] = distrib_(gen_) > cross_rate
                                ? start_population[x1_idx].weight[i]
                                : new_des.weight[i];
        if (new_des.weight[i] < new_des.min_limits[i]) {
          new_des.weight[i] = new_des.min_limits[i];
        }
        if (new_des.weight[i] > new_des.max_limits[i]) {
          new_des.weight[i] = new_des.max_limits[i];
        }
      }

      fun(&new_des);
      population.push_back(new_des);
      counter++;
    }

    start_population = population;
    output.emplace_back(*std::ranges::min_element(
        start_population,
        [&](const auto &a, const auto &b) { return a.cost < b.cost; }));


  }

  return output;
}