/*
 * Raccess-backed energy model for LinearCapR-style DP.
 */
#pragma once

#include "energy_api.hpp"
#include "../raccess/src/raccess/energy_model_api.hpp"
#include "../raccess/src/util/util.hpp"

namespace lcr {

class RaccessEnergyModel : public EnergyApi {
public:
  RaccessEnergyModel();
  void set_sequence(const std::string& seq, const std::vector<int>& seq_int) override;
  double kT() const override;
  Float energy_hairpin(int i, int j) const override;
  Float energy_loop(int i, int j, int p, int q) const override;
  Float energy_external(int i, int j) const override;
  Float energy_external_unpaired(int i, int j) const override;
  Float energy_multi_unpaired(int i, int j) const override;
  Float energy_multi_closing(int i, int j) const override;
  Float energy_multi_bif(int i, int j) const override;

private:
  typedef Raccess::ScoreModelEnergy SM;
  SM _sm;
  Raccess::EnergyModelApi _api;
};

} // namespace lcr
