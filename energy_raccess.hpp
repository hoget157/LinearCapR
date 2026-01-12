/*
 * Raccess-backed energy model for LinearCapR-style DP.
 */
#pragma once

#include "energy_api.hpp"
#include <memory>

namespace lcr {

std::unique_ptr<EnergyApi> make_raccess_energy();
std::vector<double> compute_raccess_unpaired_1(const std::string& seq, int max_span);
void debug_raccess_local(const std::string& seq, int i, int j, bool has_loop, int p, int q);

} // namespace lcr
