#pragma once

#include "types.hpp"
#include "miscs.hpp"

#include <cmath>
#include <string>

class LinearCapR{
public:
	LinearCapR(int beam_size) : beam_size(beam_size){}
	void run(const string&);
	void output(ofstream&, const string&) const;
	void clear();
// private:
	const int beam_size;

	// integerized sequence
	int seq_n;
	vector<int> seq_int;

	// DP tables
	vector<Float> alpha_O, beta_O;
	Table alpha_S, alpha_SE, alpha_M, alpha_MB, alpha_M1, alpha_M2;
	Table beta_S, beta_SE, beta_M, beta_MB, beta_M1, beta_M2;

	// calculated structural profiles
	vector<Float> prob_B, prob_I, prob_H, prob_M, prob_E, prob_S;

	// functions for pruning
	int quickselect_partition(vector<Float>&, const int, const int) const;
	Float quickselect(vector<Float>&, const int, const int, const int) const;
	Float prune(unordered_map<int, Float>&) const;

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

	inline int base_to_num(const char base) const{
		if(base == 'A' || base == 'a') return 1;
		if(base == 'C' || base == 'c') return 2;
		if(base == 'G' || base == 'g') return 3;
		if(base == 'T' || base == 't' || base == 'U' || base == 'u') return 4;
		return 0;
	}

	inline Float logsumexp(const Float x, const Float y) const{
		if(x == -INF) return y;
		if(y == -INF) return x;
		return (x > y ? x + log1p(exp(y - x)) : y + log1p(exp(x - y)));
	}

	inline Float logsumexp_equal(Float &x, const Float y) const{
		return x = logsumexp(x, y);
	}
};