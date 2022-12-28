#pragma once

#include "types.hpp"
#include "miscs.hpp"
#include "utils.hpp"

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

	// DP tables: log of sum of Boltzmann factors in interval [i, j]
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


	// returns whether base (i, j) can form pair 
	inline bool can_pair(const int i, const int j){
		return (BP_pair[seq_int[i]][seq_int[j]] > 0);
	}
	

	inline void DUMP_TABLES() const{
		DUMP("O", alpha_O);
		DUMP("S", alpha_S);
		DUMP("SE", alpha_SE);
		DUMP("M", alpha_M);
		DUMP("MB", alpha_MB);
		DUMP("M1", alpha_M1);
		DUMP("M2", alpha_M2);
	}
};