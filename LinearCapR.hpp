#pragma once

#include "types.hpp"

#include <string>

class LinearCapR{
public:
	LinearCapR(int beam_size) : beam_size(beam_size){}
	void calc_profile(const string&);
// private:
	const int beam_size;

	// integerized sequence
	vector<int> seq_int;

	// DP tables
	vector<Float> alpha_O, beta_O;
	Table alpha_S, alpha_SE, alpha_M, alpha_MB, alpha_M1, alpha_M2;
	Table beta_S, beta_SE, beta_M, beta_MB, beta_M1, beta_M2;

	int quickselect_partition(vector<Float>&, const int, const int) const;
	Float quickselect(vector<Float>&, const int, const int, const int) const;
	Float prune(unordered_map<int, Float>&);
};