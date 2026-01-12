#pragma once

#include "miscs.hpp"
#include "energy_model.hpp"

#include <string>

class LinCapR{
public:
	LinCapR(int beam_size, energy::Model model = energy::Model::Turner2004);
	void run(const string&);
	void output(ofstream&, const string&) const;
	void clear();
	Float get_energy_ensemble() const;
private:
	const energy::Params &params;
	const int beam_size;
	
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

	Float prune(Map<int, Float>&) const;

	// executable functions
	void initialize(const string &s);
	void calc_inside();
	void calc_outside();
	void calc_profile();

	// calc each energy
	Float energy_hairpin(const int, const int) const;
	Float energy_loop(const int, const int, const int, const int) const;
	Float energy_external(const int, const int) const;
	Float energy_external_unpaired(const int, const int) const;
	Float energy_multi_unpaired(const int, const int) const;
	Float energy_multi_closing(const int, const int) const;
	Float energy_multi_bif(const int, const int) const;

	int special_hairpin(const int, const int) const;

	// returns whether base (i, j) can form pair 
	inline bool can_pair(const int i, const int j) const{
		return (BP_pair[seq_int[i]][seq_int[j]] > 0);
	}
};
