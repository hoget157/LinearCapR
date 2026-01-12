/*
 * DP table utilities (re-exported) for reuse in LinearRaccess.
 */
#pragma once

#include "miscs.hpp"

namespace lcr {
namespace dp {

using ::Float;
using ::Map;
using ::Table;

inline void set_logsumexp_fast_mode() { ::set_logsumexp_fast_mode(); }
inline void set_logsumexp_legacy_mode() { ::set_logsumexp_legacy_mode(); }
inline Float logsumexp(Float x, Float y) { return ::logsumexp(x, y); }

inline Float update_sum(Table& t, const int i, const int j, const Float score) {
  return ::update_sum(t, i, j, score);
}
inline Float update_sum(vector<Float>& v, const int i, const Float score) {
  return ::update_sum(v, i, score);
}
inline Float get_value(Table& t, const int i, const int j, const Float default_value = -INF) {
  return ::get_value(t, i, j, default_value);
}
inline bool contains(const Table& t, const int i, const int j) {
  return ::contains(t, i, j);
}
inline void add_range(vector<Float>& v, const int i, const int j, const Float x) {
  ::add_range(v, i, j, x);
}
inline void prefix_sum(vector<Float>& v) { ::prefix_sum(v); }

} // namespace dp
} // namespace lcr
