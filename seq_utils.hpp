/*
 * Sequence encoding and pairing utilities (LinearCapR style).
 */
#pragma once

#include "miscs.hpp"

namespace lcr {
namespace seq {

inline int base_to_num(const char base) { return ::base_to_num(base); }

inline bool can_pair(const int a, const int b) { return (BP_pair[a][b] > 0); }

inline void build_next_pair(const vector<int>& seq_int, vector<int> next_pair[NBASE]) {
  const int seq_n = (int)seq_int.size();
  for (int i = 0; i < NBASE; ++i) next_pair[i].assign(seq_n + 1, seq_n);
  for (int i = seq_n - 1; i >= 0; --i) {
    for (int j = 0; j < NBASE; ++j) {
      next_pair[j][i] = next_pair[j][i + 1];
      if (BP_pair[seq_int[i]][j] > 0) next_pair[j][i] = i;
    }
  }
}

} // namespace seq
} // namespace lcr
