#include "energy_linearcapr.hpp"

#include <cmath>
#include <cstring>

namespace lcr {

LinearCapREnergyModel::LinearCapREnergyModel(const energy::Params& params)
	: _params(params), _seq_n(0) {}

void LinearCapREnergyModel::set_sequence(const std::string& seq, const std::vector<int>& seq_int) {
	_seq = seq;
	_seq_int = seq_int;
	_seq_n = (int)seq_int.size();
}

double LinearCapREnergyModel::kT() const { return _params.kT; }

int LinearCapREnergyModel::special_hairpin(int i, int j) const {
	if (!_params.has_special_hairpins) return -1;

	const int d = j - i - 1;
	const char *loops_seq = nullptr;
	if (d == 3) loops_seq = _params.Triloops;
	else if (d == 4) loops_seq = _params.Tetraloops;
	else if (d == 6) loops_seq = _params.Hexaloops;
	else return -1;

	if (!loops_seq) return -1;

	const std::string loop = _seq.substr(i, d + 2);
	const char *sp = strstr(loops_seq, loop.c_str());
	return (sp ? (sp - loops_seq) / (d + 3) : -1);
}

Float LinearCapREnergyModel::energy_hairpin(int i, int j) const {
	const int type = BP_pair[_seq_int[i]][_seq_int[j]];
	const int d = j - i - 1;

	const int index = special_hairpin(i, j);
	if (index != -1) {
		if (d == 3 && _params.Triloop37) return _params.Triloop37[index];
		if (d == 4 && _params.Tetraloop37) return _params.Tetraloop37[index];
		if (d == 6 && _params.Hexaloop37) return _params.Hexaloop37[index];
	}

	Float energy = (d <= MAXLOOP ? _params.hairpin37[d] : _params.hairpin37[30] + _params.lxc37 * log(d / 30.));

	if (d != 3) {
		energy += _params.mismatchH37[type][_seq_int[i + 1]][_seq_int[j - 1]];
	} else if (type > 2) {
		energy += _params.TerminalAU37;
	}
	return energy;
}

Float LinearCapREnergyModel::energy_loop(int i, int j, int p, int q) const {
	const int type1 = BP_pair[_seq_int[i]][_seq_int[j]];
	const int type2 = BP_pair[_seq_int[q]][_seq_int[p]];
	const int d1 = p - i - 1, d2 = j - q - 1;
	const int d = d1 + d2, dmin = min(d1, d2), dmax = max(d1, d2);
	const int si = _seq_int[i + 1];
	const int sj = _seq_int[j - 1];
	const int sp = _seq_int[p - 1];
	const int sq = _seq_int[q + 1];

	if (dmax == 0) {
		return _params.stack37[type1][type2];
	}

	if (dmin == 0) {
		Float energy = (d <= MAXLOOP ? _params.bulge37[d] : _params.bulge37[30] + _params.lxc37 * log(d / 30.));

		if (dmax == 1) energy += _params.stack37[type1][type2];
		else {
			if (type1 > 2) energy += _params.TerminalAU37;
			if (type2 > 2) energy += _params.TerminalAU37;
		}
		return energy;
	}

	if (d1 == 1 && d2 == 1) return (*_params.int11_37)[type1][type2][si][sj];
	if (d1 == 1 && d2 == 2) return (*_params.int21_37)[type1][type2][si][sq][sj];
	if (d1 == 2 && d2 == 1) return (*_params.int21_37)[type2][type1][sq][si][sp];
	if (d1 == 2 && d2 == 2) return (*_params.int22_37)[type1][type2][si][sp][sq][sj];

	Float energy = (d <= MAXLOOP ? _params.internal_loop37[d] : _params.internal_loop37[30] + _params.lxc37 * log(d / 30.));
	energy += min(_params.MAX_NINIO, _params.ninio37 * (dmax - dmin));

	if (dmin == 1) {
		energy += _params.mismatch1nI37[type1][si][sj] + _params.mismatch1nI37[type2][sq][sp];
	} else if (dmin == 2 && dmax == 3) {
		energy += _params.mismatch23I37[type1][si][sj] + _params.mismatch23I37[type2][sq][sp];
	} else {
		energy += _params.mismatchI37[type1][si][sj] + _params.mismatchI37[type2][sq][sp];
	}

	return energy;
}

Float LinearCapREnergyModel::energy_multi_unpaired(int i, int j) const {
	return 0;
}

Float LinearCapREnergyModel::energy_multi_closing(int i, int j) const {
	return energy_multi_bif(j, i) + _params.ML_closing37;
}

Float LinearCapREnergyModel::energy_multi_bif(int i, int j) const {
	const int type = BP_pair[_seq_int[i]][_seq_int[j]];
	Float energy = _params.ML_intern37;

	const bool has_left = (i - 1) >= 0;
	const bool has_right = (j + 1) < _seq_n;

	if (_params.allow_mismatch_multi && _params.mismatchM37 && has_left && has_right) {
		energy += _params.mismatchM37[type][_seq_int[i - 1]][_seq_int[j + 1]];
	} else {
		if (has_left) energy += _params.dangle5_37[type][_seq_int[i - 1]];
		if (has_right) energy += _params.dangle3_37[type][_seq_int[j + 1]];
	}

	if (type > 2) energy += _params.TerminalAU37;

	return energy;
}

Float LinearCapREnergyModel::energy_external(int i, int j) const {
	const int type = BP_pair[_seq_int[i]][_seq_int[j]];
	Float energy = 0;

	const bool has_left = (i - 1) >= 0;
	const bool has_right = (j + 1) < _seq_n;

	if (_params.allow_mismatch_external && _params.mismatchExt37 && has_left && has_right) {
		energy += _params.mismatchExt37[type][_seq_int[i - 1]][_seq_int[j + 1]];
	} else {
		if (has_left) energy += _params.dangle5_37[type][_seq_int[i - 1]];
		if (has_right) energy += _params.dangle3_37[type][_seq_int[j + 1]];
	}

	if (type > 2) energy += _params.TerminalAU37;

	return energy;
}

Float LinearCapREnergyModel::energy_external_unpaired(int i, int j) const {
	return 0;
}

} // namespace lcr
