#pragma once

#include "miscs.hpp"
#include "energy_model.hpp"
#include "energy_api.hpp"

#include <string>
#include <memory>

class LinCapR{
public:
	enum class EnergyEngine { LinearCapR, Raccess };
	LinCapR(int beam_size,
	       energy::Model model = energy::Model::Turner2004,
	       EnergyEngine engine = EnergyEngine::LinearCapR,
	       bool normalize_profiles = true,
	       Float normalize_warn_eps = 1e-6);
	void run(const string&);
	void output(ofstream&, const string&) const;
	void clear();
	Float get_energy_ensemble() const;
	Float get_logZ() const;
	Float get_logZ_outer_adjusted() const;
	const vector<Float>& get_prob_stem() const;
	void debug_stem_pairs(int idx, int topn) const;
	void debug_pair(int i, int j) const;
	void debug_prob(int idx) const;
	void debug_hairpin(int idx, int topn) const;
	void debug_external(int idx, int topn) const;
	void debug_internal(int idx, int topn) const;
	void debug_multi_unpaired(int i, int j) const;
	void debug_multi_prob(int idx, int topn) const;
	void debug_dp_dump(const string& state, int i0, int j0, int i1, int j1) const;
	void debug_outer_dp_dump(int i0, int i1) const;
	void set_debug_hairpin_build(bool enable);
	void set_debug_se(int i, int j, int max_hits = 200);
private:
	const energy::Params &params;
	const int beam_size;
	const bool normalize_profiles;
	const Float normalize_warn_eps;
	bool debug_hairpin_build_log = false;
	bool debug_se_log = false;
	int debug_se_i = -1;
	int debug_se_j = -1;
	int debug_se_hits = 0;
	int debug_se_max_hits = 200;
	
	string seq;

	// integerized sequence
	int seq_n;
	vector<int> seq_int;

	// next_pair[i][j] := (the first index after j (including j) which can be a pair with base i, otherwise seq_n)
	vector<int> next_pair[NBASE];

	// DP tables: log of sum of Boltzmann factors in interval [i, j]
	vector<Float> alpha_O, beta_O;
	Table alpha_S, alpha_SE, alpha_M, alpha_MB, alpha_M1, alpha_M2, *alphas[NTABLES];
	Table beta_S, beta_SE, beta_M, beta_MB, beta_M1, beta_M2, *betas[NTABLES];

	// calculated structural profiles
	vector<Float> prob_B, prob_I, prob_H, prob_M, prob_E, prob_S, *probs[NPROBS];

	std::unique_ptr<lcr::EnergyApi> _energy;

	Float prune(Map<int, Float>&) const;

	// executable functions
	void initialize(const string &s);
	void calc_inside();
	void calc_outside();
	void calc_profile();

	// returns whether base (i, j) can form pair 
	inline bool can_pair(const int i, const int j) const{
		return (BP_pair[seq_int[i]][seq_int[j]] > 0);
	}
};
