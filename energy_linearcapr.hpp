/*
 * LinearCapR energy model implementation.
 */
#pragma once

#include "energy_api.hpp"
#include "energy_model.hpp"

namespace lcr {

class LinearCapREnergyModel : public EnergyApi {
public:
  explicit LinearCapREnergyModel(const energy::Params& params);
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
  int special_hairpin(int i, int j) const;

  const energy::Params& _params;
  std::string _seq;
  std::vector<int> _seq_int;
  int _seq_n;
};

} // namespace lcr
