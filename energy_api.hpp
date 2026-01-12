/*
 * Energy API for LinearCapR-style DP.
 */
#pragma once

#include "miscs.hpp"
#include <string>
#include <vector>

namespace lcr {

class EnergyApi {
public:
  virtual ~EnergyApi() = default;
  virtual void set_sequence(const std::string& seq, const std::vector<int>& seq_int) = 0;
  virtual double kT() const = 0;
  virtual Float energy_hairpin(int i, int j) const = 0;
  virtual Float energy_loop(int i, int j, int p, int q) const = 0;
  virtual Float energy_external(int i, int j) const = 0;
  virtual Float energy_external_unpaired(int i, int j) const = 0;
  virtual Float energy_multi_unpaired(int i, int j) const = 0;
  virtual Float energy_multi_closing(int i, int j) const = 0;
  virtual Float energy_multi_bif(int i, int j) const = 0;
};

} // namespace lcr
