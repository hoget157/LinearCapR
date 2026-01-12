#include "energy_raccess.hpp"

namespace lcr {

RaccessEnergyModel::RaccessEnergyModel()
	: _sm(), _api(_sm) {
	_api.initialize();
}

void RaccessEnergyModel::set_sequence(const std::string& seq, const std::vector<int>&) {
	SM::Seq codes;
	codes.resize(seq.size());
	::Alpha::str_to_ncodes(seq.begin(), seq.end(), codes.begin());
	_api.set_seq(codes);
}

double RaccessEnergyModel::kT() const {
	return Raccess::ScoreModelEnergy::RT_KCAL_MOL();
}

Float RaccessEnergyModel::energy_hairpin(int i, int j) const {
	const int a = i + 1;
	const int b = j + 1;
	return _api.score_to_energy(_api.log_boltz_hairpin_closed(a, b));
}

Float RaccessEnergyModel::energy_loop(int i, int j, int p, int q) const {
	const int a = i + 1;
	const int b = j + 1;
	const int c = p + 1;
	const int d = q + 1;
	return _api.score_to_energy(_api.log_boltz_loop_closed(a, b, c, d));
}

Float RaccessEnergyModel::energy_external(int i, int j) const {
	const int a = i + 1;
	const int b = j + 1;
	return _api.score_to_energy(_api.log_boltz_outer_branch_closed(a, b));
}

Float RaccessEnergyModel::energy_external_unpaired(int i, int j) const {
	const int len = (j - i + 1);
	return _api.score_to_energy(_api.log_boltz_outer_extend(0, len));
}

Float RaccessEnergyModel::energy_multi_unpaired(int i, int j) const {
	const int len = (j - i + 1);
	return _api.score_to_energy(_api.log_boltz_multi_extend(0, len));
}

Float RaccessEnergyModel::energy_multi_closing(int i, int j) const {
	const int a = i + 1;
	const int b = j + 1;
	return _api.score_to_energy(_api.log_boltz_multi_close_closed(a, b));
}

Float RaccessEnergyModel::energy_multi_bif(int i, int j) const {
	const int a = i + 1;
	const int b = j + 1;
	return _api.score_to_energy(_api.log_boltz_multi_open_closed(a, b));
}

} // namespace lcr
