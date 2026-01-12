/*
 * Raccess-backed energy model for LinearCapR-style DP.
 */
#pragma once

#include "energy_api.hpp"
#include <memory>

namespace lcr {

std::unique_ptr<EnergyApi> make_raccess_energy();

} // namespace lcr
