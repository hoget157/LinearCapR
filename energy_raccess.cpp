#include "energy_raccess.hpp"

// Avoid macro conflicts (e.g. MAXLOOP) with Raccess headers.
#ifdef MAXLOOP
#pragma push_macro("MAXLOOP")
#undef MAXLOOP
#define LCR_RESTORE_MAXLOOP 1
#endif
#include "../raccess/src/raccess/energy_model_api.hpp"
#include "../raccess/src/util/util.hpp"
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

} // namespace lcr
