#include "energy_raccess.hpp"

// Avoid macro conflicts (e.g. MAXLOOP) with Raccess headers.
#ifdef MAXLOOP
#pragma push_macro("MAXLOOP")
#undef MAXLOOP
#define LCR_RESTORE_MAXLOOP 1
#endif
#include "../raccess/src/raccess/energy_model_api.hpp"
#include "../raccess/src/raccess/prob_model.hpp"
#include "../raccess/src/util/util.hpp"
#include <cstdio>
#ifdef LCR_RESTORE_MAXLOOP
#pragma pop_macro("MAXLOOP")
#undef LCR_RESTORE_MAXLOOP
#endif

namespace lcr {

class RaccessEnergyModelImpl : public EnergyApi {
public:
	RaccessEnergyModelImpl() : _sm(), _api(_sm) { _api.initialize(); }
	void set_sequence(const std::string& seq, const std::vector<int>&) override {
		SM::Seq codes;
		codes.resize(seq.size());
		::Alpha::str_to_ncodes(seq.begin(), seq.end(), codes.begin());
		_api.set_seq(codes);
	}
	double kT() const override { return Raccess::ScoreModelEnergy::RT_KCAL_MOL(); }
	Float energy_hairpin(int i, int j) const override {
		const int a = i + 1;
		const int b = j + 1;
		return _api.score_to_energy(_api.log_boltz_hairpin_closed(a, b));
	}
	Float energy_loop(int i, int j, int p, int q) const override {
		const int a = i + 1;
		const int b = j + 1;
		const int c = p + 1;
		const int d = q + 1;
		return _api.score_to_energy(_api.log_boltz_loop_closed(a, b, c, d));
	}
	Float energy_external(int i, int j) const override {
		const int a = i + 1;
		const int b = j + 1;
		return _api.score_to_energy(_api.log_boltz_outer_branch_closed(a, b));
	}
	Float energy_external_unpaired(int i, int j) const override {
		const int len = (j - i + 1);
		return _api.score_to_energy(_api.log_boltz_outer_extend(0, len));
	}
	Float energy_multi_unpaired(int i, int j) const override {
		const int len = (j - i + 1);
		return _api.score_to_energy(_api.log_boltz_multi_extend(0, len));
	}
	Float energy_multi_closing(int i, int j) const override {
		const int a = i + 1;
		const int b = j + 1;
		return _api.score_to_energy(_api.log_boltz_multi_close_closed(a, b));
	}
	Float energy_multi_bif(int i, int j) const override {
		const int a = i + 1;
		const int b = j + 1;
		return _api.score_to_energy(_api.log_boltz_multi_open_closed(a, b));
	}

private:
	typedef Raccess::ScoreModelEnergy SM;
	SM _sm;
	Raccess::EnergyModelApi _api;
};

std::unique_ptr<EnergyApi> make_raccess_energy() {
	return std::unique_ptr<EnergyApi>(new RaccessEnergyModelImpl());
}

std::vector<double> compute_raccess_unpaired_1(const std::string& seq, int max_span) {
	using SM = Raccess::ScoreModelEnergy;
	using PM = Raccess::ProbModel<SM>;
	SM sm;
	sm.initialize();

	SM::Seq codes;
	codes.resize(seq.size());
	::Alpha::str_to_ncodes(seq.begin(), seq.end(), codes.begin());
	sm.set_seq(codes);

	PM pm;
	pm.set_score_model(sm);
	pm.set_max_span(max_span);
	pm.set_prob_thr(0);
	PM::VI acc_lens;
	acc_lens.push_back(1);
	pm.set_acc_lens(acc_lens);

	std::vector<double> unpaired(seq.size(), 0.0);
	auto block = [&](int i, int w, double energy) {
		if (w != 1) return;
		const double logp = sm.energy_to_score(energy);
		unpaired[i] = exp(logp);
	};
	pm.compute_prob(block);
	return unpaired;
}

void debug_raccess_local(const std::string& seq, int i, int j, bool has_loop, int p, int q) {
	using SM = Raccess::ScoreModelEnergy;
	using Api = Raccess::EnergyModelApi;
	SM sm;
	sm.initialize();

	SM::Seq codes;
	codes.resize(seq.size());
	::Alpha::str_to_ncodes(seq.begin(), seq.end(), codes.begin());
	sm.set_seq(codes);

	Api api(sm);

	const int a = i + 1;
	const int b = j + 1;
	const int n = static_cast<int>(seq.size());

	const auto base_at = [&](int one_origin) -> char {
		if(one_origin < 1 || one_origin > n) return 'N';
		return seq[one_origin - 1];
	};

	const auto report = [&](const char* label, double via_closed, double direct) {
		const double diff = via_closed - direct;
		std::fprintf(stderr, "%s (%d,%d): closed=%g direct=%g diff=%g\n", label, i, j, via_closed, direct, diff);
	};

	report("hairpin",
	       api.score_to_energy(api.log_boltz_hairpin_closed(a, b)),
	       api.score_to_energy(api.log_boltz_hairpin(a + 1, b - 1)));

	report("multi_close",
	       api.score_to_energy(api.log_boltz_multi_close_closed(a, b)),
	       api.score_to_energy(api.log_boltz_multi_close(a + 1, b - 1)));

	report("multi_open",
	       api.score_to_energy(api.log_boltz_multi_open_closed(a, b)),
	       api.score_to_energy(api.log_boltz_multi_open(a - 1, b)));

	report("external",
	       api.score_to_energy(api.log_boltz_outer_branch_closed(a, b)),
	       api.score_to_energy(api.log_boltz_outer_branch(a - 1, b)));

	if (has_loop) {
		const int c = p + 1;
		const int d = q + 1;
		const int dp_i = a + 1;
		const int dp_j = b - 1;
		std::fprintf(stderr, "loop mapping closed(a,b,c,d)=(%d,%d,%d,%d) dp(i,j,ip,jp)=(%d,%d,%d,%d)\n",
		             a, b, c, d, dp_i, dp_j, c, d);
		std::fprintf(stderr, "  closing pair (a,b) bases=%c,%c inner pair (c,d) bases=%c,%c\n",
		             base_at(a), base_at(b), base_at(c), base_at(d));
		std::fprintf(stderr,
		             "  padded bases: seq(i)=%c seq(i+1)=%c seq(j)=%c seq(j+1)=%c seq(ip)=%c seq(ip+1)=%c seq(jp)=%c seq(jp+1)=%c\n",
		             base_at(dp_i),
		             base_at(dp_i + 1),
		             base_at(dp_j),
		             base_at(dp_j + 1),
		             base_at(c),
		             base_at(c + 1),
		             base_at(d),
		             base_at(d + 1));
		const int li = c - a;
		const int lj = b - d;
		const double interior_score = sm.score_interior(dp_i, dp_j, c, d);
		const double interior_nuc_score = sm.score_interior_nuc(dp_i, dp_j, c, d);
		std::fprintf(stderr, "  interior lens (li,lj)=(%d,%d) score=%g nuc_score=%g\n",
		             li, lj, interior_score, interior_nuc_score);
		const double via_closed = api.score_to_energy(api.log_boltz_loop_closed(a, b, c, d));
		double direct = 0.0;
		if ((c == (a + 1)) && (d == (b - 1))) {
			direct = api.score_to_energy(api.log_boltz_stack(a, b + 1));
		} else {
			direct = api.score_to_energy(api.log_boltz_interior(a, b - 1, c - 1, d));
		}
		std::fprintf(stderr, "loop (%d,%d,%d,%d): closed=%g direct=%g diff=%g\n",
		             i, j, p, q, via_closed, direct, via_closed - direct);

		if ((c == (a + 1)) && (d == (b - 1))) {
			const double stack_ab1 = api.score_to_energy(api.log_boltz_stack(a, b + 1));
			const double stack_a1b = api.score_to_energy(api.log_boltz_stack(a - 1, b));
			std::fprintf(stderr, "stack variants (%d,%d): a,b+1=%g a-1,b=%g diff=%g\n",
			             i, j, stack_ab1, stack_a1b, stack_ab1 - stack_a1b);
		}
	}
}

double compute_raccess_logz(const std::string& seq, int max_span) {
	using SM = Raccess::ScoreModelEnergy;
	using PM = Raccess::ProbModel<SM>;
	SM sm;
	sm.initialize();

	SM::Seq codes;
	codes.resize(seq.size());
	::Alpha::str_to_ncodes(seq.begin(), seq.end(), codes.begin());
	sm.set_seq(codes);

	PM pm;
	pm.set_score_model(sm);
	pm.set_max_span(max_span);
	pm.set_prob_thr(0);
	PM::VI acc_lens;
	acc_lens.push_back(1);
	pm.set_acc_lens(acc_lens);

	auto no_op = [](int, int, double) {};
	pm.compute_prob(no_op);
	return pm.partition_coeff();
}

} // namespace lcr
