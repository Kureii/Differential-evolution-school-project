#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "diff_evo.h"
#include "matplotlibcpp.h"
#include "test_functions.h"

namespace plt = matplotlibcpp;

#define F (0.5)
#define CR (0.8)
#define JCR (0.9)
#define TAU (0.1)

void WriteMarkdownHeader(const std::string &filename,
                         const std::string &func_name) {
  std::ofstream file(filename, std::ios::app);  // Append mód
  if (!file.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return;
  }

  file << "## " << func_name << "\n\n";
  file << "![](" << func_name << "_comparison.png)\n\n";
  file << "| Run | Algorithm | Best Cost | Worst Cost | Best Weights | Worst "
          "Weights | Final Cost | Final Weights |\n";
  file << "|-----|-----------|-----------|------------|--------------|---------"
          "------|------------|---------------|\n";

  file.close();
}

void WriteMarkdownSummary(const std::string &filename, int run,
                          const std::string &algo_name,
                          const std::vector<de_t> &results) {
  std::ofstream file(filename, std::ios::app);
  if (!file.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return;
  }

  auto best_it = std::ranges::min_element(
      results, [](const de_t &a, const de_t &b) { return a.cost < b.cost; });
  auto worst_it = std::ranges::max_element(
      results, [](const de_t &a, const de_t &b) { return a.cost < b.cost; });

  file << "| " << run + 1 << " | " << algo_name << " | " << best_it->cost
       << " | " << worst_it->cost << " | ";
  for (const auto &weight : best_it->weight) {
    file << weight << " ";
  }
  file << "| ";
  for (const auto &weight : worst_it->weight) {
    file << weight << " ";
  }
  file << "| " << results.back().cost << " | ";
  for (const auto &weight : results.back().weight) {
    file << weight << " ";
  }
  file << "|\n";

  file.close();
}

void PlotAverageResults(const std::string &label,
                        const std::vector<std::vector<de_t>> &all_results) {
  std::vector<double> x;
  std::vector<double> avg_y(all_results[0].size(), 0.0);

  for (size_t run = 0; run < all_results.size(); ++run) {
    for (size_t i = 0; i < all_results[run].size(); ++i) {
      avg_y[i] += all_results[run][i].cost;
    }
  }

  for (double &value : avg_y) {
    value /= all_results.size();
  }

  for (size_t i = 0; i < avg_y.size(); ++i) {
    x.push_back(i);
  }

  plt::named_plot(label, x, avg_y);
}

int main() {
  auto diff_evo = std::make_unique<DiffEvo>();
  diff_evo->GenerateInitPopulation(3, 200);
  diff_evo->AddMinLimits(-5.0);
  diff_evo->AddMaxLimits(5.0);

  std::ofstream md_file("results.md");
  if (md_file.is_open()) {
    md_file << "# Summary of Results\n\n";
    md_file.close();
  }

  constexpr int RUNS = 10;

  std::vector<std::vector<de_t>> all_schefel_rand1(RUNS);
  std::vector<std::vector<de_t>> all_sphare_rand1(RUNS);
  std::vector<std::vector<de_t>> all_rosenbrock_rand1(RUNS);
  std::vector<std::vector<de_t>> all_schefel_best1(RUNS);
  std::vector<std::vector<de_t>> all_sphare_best1(RUNS);
  std::vector<std::vector<de_t>> all_rosenbrock_best1(RUNS);
  std::vector<std::vector<de_t>> all_schefel_jde(RUNS);
  std::vector<std::vector<de_t>> all_sphare_jde(RUNS);
  std::vector<std::vector<de_t>> all_rosenbrock_jde(RUNS);

  for (int run = 0; run < RUNS; ++run) {
    // Schefel Function
    all_schefel_rand1[run] =
        diff_evo->Rand1(10000, F, CR, test_functions::Schefel);
    all_schefel_best1[run] =
        diff_evo->Best1(10000, F, CR, test_functions::Schefel);
    all_schefel_jde[run] =
        diff_evo->jDE(10000, F, JCR, TAU, TAU, test_functions::Schefel);

    // Sphere Function
    all_sphare_rand1[run] =
        diff_evo->Rand1(10000, F, CR, test_functions::Sphare);
    all_sphare_best1[run] =
        diff_evo->Best1(10000, F, CR, test_functions::Sphare);
    all_sphare_jde[run] =
        diff_evo->jDE(10000, F, JCR, TAU, TAU, test_functions::Sphare);

    // Rosenbrock Function
    all_rosenbrock_rand1[run] =
        diff_evo->Rand1(10000, F, CR, test_functions::Rosenbrock);
    all_rosenbrock_best1[run] =
        diff_evo->Best1(10000, F, CR, test_functions::Rosenbrock);
    all_rosenbrock_jde[run] =
        diff_evo->jDE(10000, F, JCR, TAU, TAU, test_functions::Rosenbrock);
  }

  WriteMarkdownHeader("results.md", "Schefel_Function");
  for (int run = 0; run < RUNS; ++run) {
    WriteMarkdownSummary("results.md", run, "Rand/1", all_schefel_rand1[run]);
  }
  for (int run = 0; run < RUNS; ++run) {
    WriteMarkdownSummary("results.md", run, "Best/1", all_schefel_best1[run]);
  }
  for (int run = 0; run < RUNS; ++run) {
    WriteMarkdownSummary("results.md", run, "jDE", all_schefel_jde[run]);
  }

  WriteMarkdownHeader("results.md", "Sphere_Function");
  for (int run = 0; run < RUNS; ++run) {
    WriteMarkdownSummary("results.md", run, "Rand/1", all_sphare_rand1[run]);
  }
  for (int run = 0; run < RUNS; ++run) {
    WriteMarkdownSummary("results.md", run, "Best/1", all_sphare_best1[run]);
  }
  for (int run = 0; run < RUNS; ++run) {
    WriteMarkdownSummary("results.md", run, "jDE", all_sphare_jde[run]);
  }

  WriteMarkdownHeader("results.md", "Rosenbrock_Function");
  for (int run = 0; run < RUNS; ++run) {
    WriteMarkdownSummary("results.md", run, "Rand/1",
                         all_rosenbrock_rand1[run]);
  }
  for (int run = 0; run < RUNS; ++run) {
    WriteMarkdownSummary("results.md", run, "Best/1",
                         all_rosenbrock_best1[run]);
  }
  for (int run = 0; run < RUNS; ++run) {
    WriteMarkdownSummary("results.md", run, "jDE", all_rosenbrock_jde[run]);
  }

  // Schefel Function
  plt::figure();
  PlotAverageResults("Rand/1", all_schefel_rand1);
  PlotAverageResults("Best/1", all_schefel_best1);
  PlotAverageResults("jDE", all_schefel_jde);
  plt::title("Schefel Function - Comparison of Algorithms");
  plt::xlabel("Iteration");
  plt::ylabel("Cost");
  plt::legend();
  plt::save("Schefel_Function_comparison.png");
  plt::clf();

  // Sphere Function
  plt::figure();
  PlotAverageResults("Rand/1", all_sphare_rand1);
  PlotAverageResults("Best/1", all_sphare_best1);
  PlotAverageResults("jDE", all_sphare_jde);
  plt::title("Sphere Function - Comparison of Algorithms");
  plt::xlabel("Iteration");
  plt::ylabel("Cost");
  plt::legend();
  plt::save("Sphere_Function_comparison.png");
  plt::clf();

  // Rosenbrock Function
  plt::figure();
  PlotAverageResults("Rand/1", all_rosenbrock_rand1);
  PlotAverageResults("Best/1", all_rosenbrock_best1);
  PlotAverageResults("jDE", all_rosenbrock_jde);
  plt::title("Rosenbrock Function - Comparison of Algorithms");
  plt::xlabel("Iteration");
  plt::ylabel("Cost");
  plt::legend();
  plt::save("Rosenbrock_Function_comparison.png");

  return 0;
}
