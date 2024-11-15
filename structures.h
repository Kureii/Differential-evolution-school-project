/*
 * Created by kureii on 11/15/24.
*/

#pragma once

#include <vector>

struct dif_evo_structure {
  double cost;
  std::vector<double> weight;
  std::vector<double> max_limits;
  std::vector<double> min_limits;
};

using de_t = dif_evo_structure;
