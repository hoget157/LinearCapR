#include "LinearCapR.hpp"
#include "utils.hpp"
#include "energy_param.hpp"
#include "intloops.hpp"

#include <fstream>

// partite [lower, upper) and small scores are in [lower, split)
int LinearCapR::quickselect_partition(vector<Float> &scores, const int lower, const int upper) const{
	Float pivot = scores[upper - 1];
	int i = lower, j = upper - 1;
	while(i < j){
		while (scores[i] < pivot) i++;
        while (scores[j] > pivot) j--;
        if (scores[i] == scores[j]) i++;
        else if (i < j) swap(scores[i], scores[j]);
	}
	return j;
}


// returns k-th(0-indexed) smalled score in [lower, upper)
Float LinearCapR::quickselect(vector<Float> &scores, const int lower, const int upper, const int k) const{
	DUMP(scores, lower, upper, k);
	if(upper - lower == 1) return scores[lower];
	int split = quickselect_partition(scores, lower, upper);
	int length = split - lower + 1;
	if(length == k) return scores[split];
	if(k < length) return quickselect(scores, lower, split, k);
	return quickselect(scores, split + 1, upper, k - length);
}


// prune top-k states
Float LinearCapR::prune(unordered_map<int, Float> &states) const{
	// extract scores
	vector<Float> scores;
	for(auto [i, score] : states){
		// bias
		Float new_score = (i >= 1 ? alpha_O[i - 1] : Float(0)) + score;
        scores.push_back(new_score);
	}

	// threshold
	Float threshold = quickselect(scores, 0, scores.size(), scores.size() - beam_size);

	// erase low-scored states
	for(auto it = states.begin(); it != states.end();){
		if(it->second < threshold) it = states.erase(it);
		else it++;
	}

	return threshold;
}


// output structural profile
void LinearCapR::output(ofstream &ofs, const string &seq_name) const{
	ofs << ">" + seq_name << endl;

	ofs << "Bulge ";
	for(int i = 0; i < seq_n; i++) ofs << prob_B[i] << " ";
	ofs << endl;

	ofs << "Exterior ";
	for(int i = 0; i < seq_n; i++) ofs << prob_E[i] << " ";
	ofs << endl;

	ofs << "Hairpin ";
	for(int i = 0; i < seq_n; i++) ofs << prob_H[i] << " ";
	ofs << endl;

	ofs << "Internal ";
	for(int i = 0; i < seq_n; i++) ofs << prob_I[i] << " ";
	ofs << endl;

	ofs << "Multibranch ";
	for(int i = 0; i < seq_n; i++) ofs << prob_M[i] << " ";
	ofs << endl;

	ofs << "Stem ";
	for(int i = 0; i < seq_n; i++) ofs << prob_S[i] << " ";
	ofs << endl;

	ofs << endl;
}


// clear temp tables & profiles
void LinearCapR::clear(){
	seq_int.clear();
	seq_n = 0;

	alpha_O.clear();
	alpha_S.clear();
	alpha_SE.clear();
	alpha_M.clear();
	alpha_MB.clear();
	alpha_M1.clear();
	alpha_M2.clear();
	
	beta_O.clear();
	beta_S.clear();
	beta_SE.clear();
	beta_M.clear();
	beta_MB.clear();
	beta_M1.clear();
	beta_M2.clear();

	prob_B.clear();
	prob_I.clear();
	prob_H.clear();
	prob_M.clear();
	prob_E.clear();
	prob_S.clear();
}


// calc structural profile
void LinearCapR::run(const string &seq){
	initialize(seq);
	calc_inside();
	calc_outside();
	calc_profile();
}


// initialize
void LinearCapR::initialize(const string &seq){
	// integerize sequence
	seq_n = seq.length();
	seq_int.resize(seq_n + 1);
	for(int i = 0; i < seq_n; i++){
		seq_int[i + 1] = base_to_num(seq[i]);
	}

	// prepare DP tables
	alpha_O.resize(seq_n + 1);
	alpha_S.resize(seq_n + 1);
	alpha_SE.resize(seq_n + 1);
	alpha_M.resize(seq_n + 1);
	alpha_MB.resize(seq_n + 1);
	alpha_M1.resize(seq_n + 1);
	alpha_M2.resize(seq_n + 1);

	beta_O.resize(seq_n + 1);
	beta_S.resize(seq_n + 1);
	beta_SE.resize(seq_n + 1);
	beta_M.resize(seq_n + 1);
	beta_MB.resize(seq_n + 1);
	beta_M1.resize(seq_n + 1);
	beta_M2.resize(seq_n + 1);

	prob_B.resize(seq_n);
	prob_E.resize(seq_n);
	prob_H.resize(seq_n);
	prob_I.resize(seq_n);
	prob_M.resize(seq_n);
	prob_S.resize(seq_n);
}


// calc inside variables
void LinearCapR::calc_inside(){

}


// calc outside variables
void LinearCapR::calc_outside(){

}


// calc structural profile
void LinearCapR::calc_profile(){

}


// calc energy of hairpin loop [i, j]
Float LinearCapR::energy_hairpin(int i, int j) const{
	// int type = BP_pair[seq_int[i]][seq_int[j]];
	// int d = j - i - 1;

	// // initiation
	// double energy = (d <= MAXLOOP ? hairpin37[]);
	return 0;
}