/*
 * Beam pruning utilities extracted from LinCapR.
 */
#pragma once

#include "miscs.hpp"
#include <vector>

namespace lcr {
namespace beam {

inline int quickselect_partition(vector<Float>& scores, const int lower, const int upper) {
  const Float pivot = scores[upper - 1];
  int i = lower, j = upper - 1;
  while (i < j) {
    while (scores[i] < pivot) i++;
    while (scores[j] > pivot) j--;
    if (scores[i] == scores[j]) i++;
    else if (i < j) swap(scores[i], scores[j]);
  }
  return j;
}

inline Float quickselect(vector<Float>& scores, const int lower, const int upper, const int k) {
  if (upper - lower == 1) return scores[lower];
  const int split = quickselect_partition(scores, lower, upper);
  const int length = split - lower + 1;
  if (length == k) return scores[split];
  if (k < length) return quickselect(scores, lower, split, k);
  return quickselect(scores, split + 1, upper, k - length);
}

template <typename BiasFn>
inline Float prune_states(Map<int, Float>& states, const int beam_size, BiasFn bias) {
  if (beam_size == 0 || (int)states.size() <= beam_size) return -INF;

  vector<Float> scores;
  scores.reserve(states.size());
  for (const auto [i, score] : states) {
    scores.push_back(bias(i, score));
  }

  const Float threshold = quickselect(scores, 0, scores.size(), scores.size() - beam_size);

  for (auto it = states.begin(); it != states.end();) {
    const auto [i, score] = *it;
    const Float new_score = bias(i, score);
    if (new_score <= threshold) {
      it = states.erase(it);
    } else {
      ++it;
    }
  }
  return threshold;
}

} // namespace beam
} // namespace lcr
