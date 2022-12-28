#include "LinearCapR.hpp"
#include "utils.hpp"
#include "energy_param.hpp"
#include "intloops.hpp"

#include <fstream>
#include <cmath>

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
	if(upper - lower == 1) return scores[lower];
	int split = quickselect_partition(scores, lower, upper);
	int length = split - lower + 1;
	if(length == k) return scores[split];
	if(k < length) return quickselect(scores, lower, split, k);
	return quickselect(scores, split + 1, upper, k - length);
}


// prune top-k states
Float LinearCapR::prune(unordered_map<int, Float> &states) const{
	if(beam_size == 0 || (int)states.size() <= beam_size) return -INF;
	
	// extract scores
	vector<Float> scores;
	for(const auto [i, score] : states){
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
	seq_int.resize(seq_n);
	for(int i = 0; i < seq_n; i++){
		seq_int[i] = base_to_num(seq[i]);
	}

	// prepare DP tables
	alpha_O.resize(seq_n, -INF);
	alpha_S.resize(seq_n);
	alpha_SE.resize(seq_n);
	alpha_M.resize(seq_n);
	alpha_MB.resize(seq_n);
	alpha_M1.resize(seq_n);
	alpha_M2.resize(seq_n);

	beta_O.resize(seq_n, -INF);
	beta_S.resize(seq_n);
	beta_SE.resize(seq_n);
	beta_M.resize(seq_n);
	beta_MB.resize(seq_n);
	beta_M1.resize(seq_n);
	beta_M2.resize(seq_n);

	prob_B.resize(seq_n);
	prob_E.resize(seq_n);
	prob_H.resize(seq_n);
	prob_I.resize(seq_n);
	prob_M.resize(seq_n);
	prob_S.resize(seq_n);
}


// calc inside variables
void LinearCapR::calc_inside(){
	alpha_O[0] = 0;

	for(int j = 0; j < seq_n; j++){
		// S
		prune(alpha_S[j]);
		for(const auto &[i, score] : alpha_S[j]){
			// S -> S
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				update_sum(alpha_S, i - 1, j + 1, score - energy_loop(i - 1, j + 1, i, j) / kT);
			}
			
			// M2 -> S
			for(int n = 0; n <= MAXLOOP; n++){
				if(j + n >= seq_n) continue;
				update_sum(alpha_M2, i, j + n, score - (energy_multi_bif(i, j) + energy_multi_unpaired(j + 1, j + n)) / kT);
			}

			// SE -> S: p..i..j..q
			for(int p = i; i - p <= MAXLOOP; p--){
				for(int q = j; (q - j) + (i - p) <= MAXLOOP; q++){
					if((p == i && q == j) || p - 1 < 0 || q + 1 >= seq_n || !can_pair(p - 1, q + 1)) continue;
					update_sum(alpha_SE, p, q, score - energy_loop(p - 1, q + 1, i, j) / kT);
				}
			}

			// O -> O + S
			update_sum(alpha_O, j, (i - 1 >= 0 ? alpha_O[i - 1] : 0) + score - energy_external(i, j) / kT);
		}

		// M2
		prune(alpha_M2[j]);
		for(const auto &[i, score] : alpha_M2[j]){
			// M1 -> M2
			update_sum(alpha_M1, i, j, score);

			// MB -> M1 + M2
			if(i - 1 >= 0){
				for(const auto &[k, score_m1] : alpha_M1[i - 1]){
					update_sum(alpha_MB, k, j, score_m1 + score);
				}
			}
		}

		// MB
		prune(alpha_MB[j]);
		for(const auto &[i, score] : alpha_MB[j]){
			// M1 -> MB
			update_sum(alpha_M1, i, j, score);

			// M -> MB
			for(int n = 0; n <= MAXLOOP; n++){
				if(i - n < 0) continue;
				update_sum(alpha_M, i - n, j, score);
			}
		}

		// M1
		prune(alpha_M1[j]);

		// M
		prune(alpha_M[j]);
		for(const auto &[i, score] : alpha_M[j]){
			// SE -> M
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				update_sum(alpha_SE, i, j, score - energy_multi_closing(i - 1, j + 1) / kT);
			}
		}

		// SE -> (Hairpin)
		for(int n = TURN; n <= MAXLOOP; n++){
			int i = j - n + 1;
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				update_sum(alpha_SE, i, j, -energy_hairpin(i - 1, j + 1) / kT);
			}
		}

		// SE
		prune(alpha_SE[j]);
		for(const auto &[i, score] : alpha_SE[j]){
			// S -> SE
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				update_sum(alpha_S, i - 1, j + 1, score);
			}
		}

		// O
		if(j + 1 < seq_n){
			update_sum(alpha_O, j + 1, alpha_O[j] - energy_external_unpaired(j + 1, j + 1) / kT);
		}
	}
	DUMP_TABLES();
}


// calc outside variables
void LinearCapR::calc_outside(){

}


// calc structural profile
void LinearCapR::calc_profile(){

}


// calc energy of hairpin loop [i, j]
Float LinearCapR::energy_hairpin(int i, int j) const{
	int type = BP_pair[seq_int[i]][seq_int[j]];
	int d = j - i - 1;

	// initiation
	Float energy = (d <= MAXLOOP ? hairpin37[d] : hairpin37[30] + lxc37 * log(d / 30.));
	
	if(d != 3){
		energy += mismatchH37[type][seq_int[i + 1]][seq_int[j - 1]];
	}else if(type > 2){
		// TODO: why add terminal penalty?
		energy += TerminalAU37;
	}
	return energy;
}


// calc energy of loop [i, p, q, j]
Float LinearCapR::energy_loop(int i, int j, int p, int q) const{
	int type1 = BP_pair[seq_int[i]][seq_int[j]], type2 = BP_pair[seq_int[q]][seq_int[p]];;
	int d1 = p - i - 1, d2 = j - q - 1;
	int d = d1 + d2, dmin = min(d1, d2), dmax = max(d1, d2);
	int si = seq_int[i + 1];
	int sj = seq_int[j - 1];
	int sp = seq_int[p - 1];
	int sq = seq_int[q + 1];

	if(dmax == 0){
		// stack
		return stack37[type1][type2];
	}

	if(dmin == 0){
		// bulge
		Float energy = (d <= MAXLOOP ? bulge37[d] : bulge37[30] + lxc37 * log(d / 30.));

		if(dmax == 1) energy += stack37[type1][type2];
		else{
			if(type1 > 2) energy += TerminalAU37;
			if(type2 > 2) energy += TerminalAU37;
		}
		return energy;
	}

	// internal
	// specieal internal loops
	if(d1 == 1 && d2 == 1) return int11_37[type1][type2][si][sj];
	if(d1 == 1 && d2 == 2) return int21_37[type1][type2][si][sq][sj];
	if(d1 == 2 && d2 == 1) return int21_37[type2][type1][sq][si][sp];
	if(d1 == 2 && d2 == 2) return int22_37[type1][type2][si][sp][sq][sj];

	// generic internal loop
	Float energy =  (d <= MAXLOOP ? internal_loop37[d] : internal_loop37[30] + lxc37 * log(d / 30.));
	energy += min(MAX_NINIO, ninio37 * (dmax - dmin));
	
	// mismatch: different for sizes
	if(dmin == 1){ // 1xn
		energy += mismatch1nI37[type1][si][sj] + mismatch1nI37[type2][sq][sp];
	}else if(dmin == 2 && dmax == 3){ // 2x3
		energy += mismatch23I37[type1][si][sj] + mismatch23I37[type2][sq][sp];
	}else{ // others
		energy += mismatchI37[type1][si][sj] + mismatchI37[type2][sq][sp];
	}

	return energy;
}


// calc energy where bases in multi [i, j] are unpaired
Float LinearCapR::energy_multi_unpaired(const int i, const int j) const{
	return 0;
}


// calc energy of multiloop [i, j]
Float LinearCapR::energy_multi_closing(const int i, const int j) const{
	// we look clockwise, so i, j are swapped
	return energy_multi_bif(j, i) + ML_closing37;
}


// calc energy of bifurcation [i, j] in a multiloop
Float LinearCapR::energy_multi_bif(const int i, const int j) const{
	int type = BP_pair[seq_int[i]][seq_int[j]];
	Float energy = ML_intern37;

	if(i - 1 >= 0 && j + 1 < seq_n) energy += mismatchM37[type][seq_int[i - 1]][seq_int[j + 1]];
	else if(i - 1 >= 0) energy += dangle5_37[type][seq_int[i - 1]];
	else if(j + 1 < seq_n) energy += dangle3_37[type][seq_int[j + 1]];

	if(type > 2) energy += TerminalAU37;

	return energy;
}


// calc energy of external loop [i, j]
Float LinearCapR::energy_external(const int i, const int j) const{
	int type = BP_pair[seq_int[i]][seq_int[j]];
	Float energy = 0;

	if(i - 1 >= 0 && j + 1 < seq_n) energy += mismatchExt37[type][seq_int[i - 1]][seq_int[j + 1]];
	else if(i - 1 >= 0) energy += dangle5_37[type][seq_int[i - 1]];
	else if(j + 1 < seq_n) energy += dangle3_37[type][seq_int[j + 1]];

	if(type > 2) energy += TerminalAU37;

	return energy;
}


// calc energy where bases in external [i, j] are unpaired
Float LinearCapR::energy_external_unpaired(const int i, const int j) const{
	return 0;
}