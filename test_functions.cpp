/*
 * Created by kureii on 11/15/24.
*/

#include "test_functions.h"

#include <algorithm>

#include <cfloat>
#include <cmath>

namespace test_functions {

void Rosenbrock(de_t *des) {
  if (des->weight.size() > 1) {
    if (des->cost != 0) des->cost = 0;
    for (size_t i = 0; i < des->weight.size() - 1; ++i) {
      const auto x = des->weight[i];
      const auto x_next = des->weight[i + 1];
      des->cost = 100 * pow(x_next - pow(x, 2), 2) + pow(x - 1, 2);
    }
  } else {
    des->cost = DBL_MAX;
  }
}

void Sphare(de_t *des) {
  if (!des->weight.empty()) {
    if (des->cost != 0) des->cost = 0;
    std::ranges::for_each(des->weight, [&des](auto &v) {
        des->cost += pow(v, 2);
    });
  } else {
    des->cost = DBL_MAX;
  }
}

void Schefel(de_t *des) {
  if (!des->weight.empty()) {
    if (des->cost != 0) des->cost = 0;
    std::ranges::for_each(des->weight, [&des](auto &v) {
        des->cost += sin(sqrt(fabs(v)));
    });
  } else {
    des->cost = DBL_MAX;
  }
}

}