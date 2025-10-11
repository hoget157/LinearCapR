#include "LinCapR.hpp"

#include <fstream>
#include <algorithm>
#include <cstring>

LinCapR::LinCapR(int beam_size, energy::Model model)
	: params(energy::get_params(model)), beam_size(beam_size){}


// partite [lower, upper) and small scores are in [lower, split]
int LinCapR::quickselect_partition(vector<Float> &scores, const int lower, const int upper) const{
	const Float pivot = scores[upper - 1];
	int i = lower, j = upper - 1;
	while(i < j){
		while (scores[i] < pivot) i++;
        while (scores[j] > pivot) j--;
        if (scores[i] == scores[j]) i++;
        else if (i < j) swap(scores[i], scores[j]);
	}
	return j;
}


// returns k-th(1-indexed) smallest score in [lower, upper)
Float LinCapR::quickselect(vector<Float> &scores, const int lower, const int upper, const int k) const{
	if(upper - lower == 1) return scores[lower];
	const int split = quickselect_partition(scores, lower, upper);
	const int length = split - lower + 1;
	if(length == k) return scores[split];
	if(k < length) return quickselect(scores, lower, split, k);
	return quickselect(scores, split + 1, upper, k - length);
}


// prune top-k states
Float LinCapR::prune(Map<int, Float> &states) const{
	if(beam_size == 0 || (int)states.size() <= beam_size) return -INF;
	
	// extract scores
	vector<Float> scores;
	for(const auto [i, score] : states){
		// bias
		const Float new_score = (i >= 1 ? alpha_O[i - 1] : Float(0)) + score;
        scores.push_back(new_score);
	}

	// threshold
	const Float threshold = quickselect(scores, 0, scores.size(), scores.size() - beam_size);

	// erase low-scored states
	for(auto it = states.begin(); it != states.end();){
		const auto [i, score] = *it;
		const Float new_score = (i >= 1 ? alpha_O[i - 1] : Float(0)) + score;
		if(new_score <= threshold){
			// unorderd_map
			it = states.erase(it);
			// google hash
			// states.erase(it++);
		}
		else it++;
	}
	// google hash
	// states.resize(0);

	return threshold;
}


// output structural profile
void LinCapR::output(ofstream &ofs, const string &seq_name) const{
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
void LinCapR::clear(){
	seq = "";
	seq_int.clear();
	seq_n = 0;

	for(int i = 0; i < NTABLES; i++){
		alphas[i]->clear();
		betas[i]->clear();
	}

	for(int i = 0; i < NPROBS; i++) probs[i]->clear();
}


// returns free energy of ensemble in kcal/mol
Float LinCapR::get_energy_ensemble() const{
	return (alpha_O[seq_n - 1] * -(params.temperature + params.k0) * params.gas_constant) / 1000;
}


// calc structural profile
void LinCapR::run(const string &seq){
	initialize(seq);
	calc_inside();
	calc_outside();
	calc_profile();
}


// initialize
void LinCapR::initialize(const string &seq){
	this->seq = seq;

	// integerize sequence
	seq_n = seq.length();
	seq_int.resize(seq_n);
	for(int i = 0; i < seq_n; i++){
		seq_int[i] = base_to_num(seq[i]);
	}

	// prepare DP tables
	alphas[0] = &alpha_S;
	alphas[1] = &alpha_SE;
	alphas[2] = &alpha_M;
	alphas[3] = &alpha_MB;
	alphas[4] = &alpha_M1;
	alphas[5] = &alpha_M2;

	betas[0] = &beta_S;
	betas[1] = &beta_SE;
	betas[2] = &beta_M;
	betas[3] = &beta_MB;
	betas[4] = &beta_M1;
	betas[5] = &beta_M2;

	alpha_O.resize(seq_n, -INF);
	beta_O.resize(seq_n, -INF);
	for(int i = 0; i < NTABLES; i++){
		alphas[i]->resize(seq_n);
		betas[i]->resize(seq_n);
		// google hash
		// for(int j = 0; j < seq_n; j++){
		// 	alphas[i]->at(j).set_empty_key(-1);
		// 	alphas[i]->at(j).set_deleted_key(-2);
		// 	betas[i]->at(j).set_empty_key(-1);
		// 	betas[i]->at(j).set_deleted_key(-2);
		// }
	}

	// prepare prob vectors
	probs[0] = &prob_B;
	probs[1] = &prob_E;
	probs[2] = &prob_H;
	probs[3] = &prob_I;
	probs[4] = &prob_M;
	probs[5] = &prob_S;

	for(int i = 0; i < NPROBS; i++) probs[i]->resize(seq_n);

	// calc next pair index
	for(int i = 0; i < NBASE; i++) next_pair[i].resize(seq_n + 1, seq_n);
	for(int i = seq_n - 1; i >= 0; i--){
		for(int j = 0; j < NBASE; j++){
			next_pair[j][i] = next_pair[j][i + 1];
			if(BP_pair[seq_int[i]][j] > 0) next_pair[j][i] = i;
		}
	}
}


// calc inside variables
void LinCapR::calc_inside(){
	alpha_O[0] = 0;

	for(int j = 0; j < seq_n; j++){
		// S
		prune(alpha_S[j]);
		for(const auto [i, score] : alpha_S[j]){
			// S -> S
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				update_sum(alpha_S, i - 1, j + 1, score - energy_loop(i - 1, j + 1, i, j) / params.kT);
			}
			
			// M2 -> S
			for(int n = 0; n <= MULTI_MAX_UNPAIRED && j + n < seq_n; n++){
				update_sum(alpha_M2, i, j + n, score - (energy_multi_bif(i, j) + energy_multi_unpaired(j + 1, j + n)) / params.kT);
			}

			// SE -> S: p..i..j..q, [p - 1, q] can be pair
			for(int p = i; i - p <= MAXLOOP && p >= 1; p--){
				for(int q = next_pair[seq_int[p - 1]][j + 1]; q < seq_n && (q - j - 1) + (i - p) <= MAXLOOP; q = next_pair[seq_int[p - 1]][q + 1]){
					if((p == i && q == j + 1)) continue;
					update_sum(alpha_SE, p, q - 1, score - energy_loop(p - 1, q, i, j) / params.kT);
				}
			}

			// O -> O + S
			update_sum(alpha_O, j, (i - 1 >= 0 ? alpha_O[i - 1] : 0) + score - energy_external(i, j) / params.kT);
		}

		// M2
		prune(alpha_M2[j]);
		for(const auto [i, score] : alpha_M2[j]){
			// M1 -> M2
			update_sum(alpha_M1, i, j, score);

			// MB -> M1 + M2
			if(i - 1 >= 0){
				for(const auto [k, score_m1] : alpha_M1[i - 1]){
					update_sum(alpha_MB, k, j, score_m1 + score);
				}
			}
		}

		// MB
		prune(alpha_MB[j]);
		for(const auto [i, score] : alpha_MB[j]){
			// M1 -> MB
			update_sum(alpha_M1, i, j, score);

			// M -> MB
			for(int n = 0; n <= MULTI_MAX_UNPAIRED && i - n >= 0; n++){
				update_sum(alpha_M, i - n, j, score);
			}
		}

		// M1
		prune(alpha_M1[j]);

		// M
		prune(alpha_M[j]);
		for(const auto [i, score] : alpha_M[j]){
			// SE -> M
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				update_sum(alpha_SE, i, j, score - energy_multi_closing(i - 1, j + 1) / params.kT);
			}
		}

		// SE -> (Hairpin)
		for(int n = TURN; n <= MAXLOOP; n++){
			const int i = j - n + 1;
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				update_sum(alpha_SE, i, j, -energy_hairpin(i - 1, j + 1) / params.kT);
			}
		}

		// SE
		prune(alpha_SE[j]);
		for(const auto [i, score] : alpha_SE[j]){
			// S -> SE
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				update_sum(alpha_S, i - 1, j + 1, score);
			}
		}

		// O -> O
		if(j + 1 < seq_n){
			update_sum(alpha_O, j + 1, alpha_O[j] - energy_external_unpaired(j + 1, j + 1) / params.kT);
		}
	}
}


// calc outside variables
void LinCapR::calc_outside(){
	for(int j = seq_n - 1; j >= 0; j--){
		// O
		// O -> O
		update_sum(beta_O, j, (j + 1 < seq_n ? beta_O[j + 1] : 0) - energy_external_unpaired(j + 1, j + 1) / params.kT);
		
		// O -> O + S
		for(const auto [i, score] : alpha_S[j]){
			update_sum(beta_O, i, score + (j + 1 < seq_n ? beta_O[j + 1] : 0) - energy_external(i, j) / params.kT);
		}

		// SE
		for(const auto [i, _] : alpha_SE[j]){
			// S -> SE
			if(i - 1 >= 0 && j + 1 < seq_n){
				update_sum(beta_SE, i, j, get_value(beta_S, i - 1, j + 1));
			}
		}

		// M
		for(const auto [i, _] : alpha_M[j]){
			// SE -> M
			if(i - 1 >= 0 && j + 1 < seq_n){
				update_sum(beta_M, i, j, get_value(beta_SE, i, j) - energy_multi_closing(i - 1, j + 1) / params.kT);
			}
		}

		// MB
		for(const auto [i, _] : alpha_MB[j]){
			// M1 -> MB
			update_sum(beta_MB, i, j, get_value(beta_M1, i, j));

			// M -> MB
			for(int n = 0; n <= MULTI_MAX_UNPAIRED && i - n >= 0; n++){
				update_sum(beta_MB, i, j, get_value(beta_M, i - n, j));
			}
		}

		// M1, M2
		for(const auto [i, score_M2] : alpha_M2[j]){
			// M1 -> M2
			update_sum(beta_M2, i, j, get_value(beta_M1, i, j));

			// MB -> M1 + M2
			if(i - 1 < 0) continue;
			for(const auto [k, score_M1] : alpha_M1[i - 1]){
				update_sum(beta_M1, k, i - 1, get_value(beta_MB, k, j) + score_M2);
				update_sum(beta_M2, i, j, get_value(beta_MB, k, j) + score_M1);
			}
		}

		// S
		for(const auto [i, _] : alpha_S[j]){
			// O -> O + S
			update_sum(beta_S, i, j, (i - 1 >= 0 ? alpha_O[i - 1] : 0) + (j + 1 < seq_n ? beta_O[j + 1] : 0) - energy_external(i, j) / params.kT);

			// SE -> S
			for(int p = i; i - p <= MAXLOOP && p >= 1; p--){
				for(int q = next_pair[seq_int[p - 1]][j + 1]; q < seq_n && (q - j - 1) + (i - p) <= MAXLOOP; q = next_pair[seq_int[p - 1]][q + 1]){
					if((p == i && q == j + 1)) continue;
					update_sum(beta_S, i, j, get_value(beta_SE, p, q - 1) - energy_loop(p - 1, q, i, j) / params.kT);
				}
			}

			// S -> S
			if(i - 1 >= 0 && j + 1 < seq_n){
				update_sum(beta_S, i, j, get_value(beta_S, i - 1, j + 1) - energy_loop(i - 1, j + 1, i, j) / params.kT);
			}
			
			// M2 -> S
			for(int n = 0; n <= MULTI_MAX_UNPAIRED && j + n < seq_n; n++){
				update_sum(beta_S, i, j, get_value(beta_M2, i, j + n) - (energy_multi_bif(i, j) + energy_multi_unpaired(j + 1, j + n)) / params.kT);
			}
		}
	}
}


// calc structural profile
void LinCapR::calc_profile(){
	const Float logZ = alpha_O[seq_n - 1];

	for(int k = 0; k < seq_n; k++){
		for(const auto [j, score] : beta_SE[k]){
			// H
			add_range(prob_H, j, k, exp(score - energy_hairpin(j - 1, k + 1) / params.kT - logZ));

			// B, I
			for(int p = j; p <= min(j + MAXLOOP, k - 1); p++){
				for(int q = k; q >= p + TURN + 1 && (p - j) + (k - q) <= MAXLOOP; q--){
					if((p == j && q == k) || !contains(alpha_S, p, q)) continue;
					const Float new_score = exp(score + alpha_S[q][p] - energy_loop(j - 1, k + 1, p, q) / params.kT - logZ);
					add_range((q == k ? prob_B : prob_I), j, p - 1, new_score);
					add_range((p == j ? prob_B : prob_I), q + 1, k, new_score);
				}
			}
		}
	}
	prefix_sum(prob_B);
	prefix_sum(prob_H);
	prefix_sum(prob_I);

	// M
	for(int k = 0; k < seq_n; k++){
		for(const auto [p, score] : alpha_MB[k]){
			for(int j = p - 1; j >= max(0, p - MAXLOOP); j--){
				if(!contains(beta_M, j, k)) continue;
				const Float new_score = exp(score + beta_M[k][j] - energy_multi_unpaired(j, p - 1) / params.kT - logZ);
				add_range(prob_M, j, p - 1, new_score);
			}
		}
	}
	for(int q = 0; q < seq_n; q++){
		for(const auto [j, score] : alpha_S[q]){
			for(int k = q + 1; k <= min(seq_n - 1, q + MAXLOOP); k++){
				if(!contains(beta_M2, j, k)) continue;
				const Float new_score = exp(score + beta_M2[k][j] - (energy_multi_bif(j, q) + energy_multi_unpaired(q + 1, k)) / params.kT - logZ);
				add_range(prob_M, q + 1, k, new_score);
			}
		}
	}
	prefix_sum(prob_M);

	// S
	for(int j = 0; j < seq_n; j++){
		for(const auto [i, score] : alpha_S[j]){
			const Float new_score = exp(score + beta_S[j][i] - logZ);
			prob_S[i] += new_score;
			prob_S[j] += new_score;
		}
	}

	// E
	prob_E[0] = exp(beta_O[1] - logZ);
	prob_E[seq_n - 1] = exp(alpha_O[seq_n - 2] - logZ);
	for(int i = 1; i < seq_n - 1; i++){
		prob_E[i] = exp(alpha_O[i - 1] + beta_O[i + 1] - logZ);
	}

	// regularize
	for(int i = 0; i < seq_n; i++){
		Float sum_prob_i = 0;
		for(int j = 0; j < NPROBS; j++){
			// negative probabilities to 0
			if(probs[j]->at(i) < 0) probs[j]->at(i) = 0;
			sum_prob_i += probs[j]->at(i);
		}

		// sum of probabilities to 1
		for(int j = 0; j < NPROBS; j++) probs[j]->at(i) /= sum_prob_i;
	}
}


// returns index if loop [i, j] is special hairpin, otherwise -1
int LinCapR::special_hairpin(const int i, const int j) const{
#ifdef LEGACY_ENERGY
	return -1;
#else
	const int d = j - i - 1;
	const char *loops_seq;
	if(d == 3) loops_seq = Triloops;
	else if(d == 4) loops_seq = Tetraloops;
	else if(d == 6) loops_seq = Hexaloops;
	else return -1;

	const char *sp = strstr(loops_seq, seq.substr(i, d + 2).c_str());
	return (sp ? (sp - loops_seq) / (d + 3) : -1);
#endif
}


// calc energy of hairpin loop [i, j]
Float LinCapR::energy_hairpin(const int i, const int j) const{
	const int type = BP_pair[seq_int[i]][seq_int[j]];
	const int d = j - i - 1;
	
#ifndef LEGACY_ENERGY
	// check special hairpin
	const int index = special_hairpin(i, j);
	if(index != -1){
		if(d == 3) return Triloop37[index];
		if(d == 4) return Tetraloop37[index];
		if(d == 6) return Hexaloop37[index];
	}
#endif

	// initiation
	Float energy = (d <= MAXLOOP ? hairpin37[d] : hairpin37[30] + lxc37 * log(d / 30.));
	
	if(d != 3){
		energy += mismatchH37[type][seq_int[i + 1]][seq_int[j - 1]];
	}else if(type > 2){
		energy += TerminalAU37;
	}
	return energy;
}


// calc energy of loop [i, p, q, j]
Float LinCapR::energy_loop(const int i, const int j, const int p, const int q) const{
	const int type1 = BP_pair[seq_int[i]][seq_int[j]], type2 = BP_pair[seq_int[q]][seq_int[p]];;
	const int d1 = p - i - 1, d2 = j - q - 1;
	const int d = d1 + d2, dmin = min(d1, d2), dmax = max(d1, d2);
	const int si = seq_int[i + 1];
	const int sj = seq_int[j - 1];
	const int sp = seq_int[p - 1];
	const int sq = seq_int[q + 1];

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
Float LinCapR::energy_multi_unpaired(const int i, const int j) const{
	return 0;
}


// calc energy of multiloop [i, j]
Float LinCapR::energy_multi_closing(const int i, const int j) const{
	// we look clockwise, so i, j are swapped
	return energy_multi_bif(j, i) + ML_closing37;
}


// calc energy of bifurcation [i, j] in a multiloop
Float LinCapR::energy_multi_bif(const int i, const int j) const{
	const int type = BP_pair[seq_int[i]][seq_int[j]];
	Float energy = ML_intern37;

#ifdef LEGACY_ENERGY
	if(i - 1 >= 0) energy += dangle5_37[type][seq_int[i - 1]];
	if(j + 1 < seq_n) energy += dangle3_37[type][seq_int[j + 1]];
#else
	if(i - 1 >= 0 && j + 1 < seq_n) energy += mismatchM37[type][seq_int[i - 1]][seq_int[j + 1]];
	else if(i - 1 >= 0) energy += dangle5_37[type][seq_int[i - 1]];
	else if(j + 1 < seq_n) energy += dangle3_37[type][seq_int[j + 1]];
#endif

	if(type > 2) energy += TerminalAU37;

	return energy;
}


// calc energy of external loop [i, j]
Float LinCapR::energy_external(const int i, const int j) const{
	const int type = BP_pair[seq_int[i]][seq_int[j]];
	Float energy = 0;

#ifdef LEGACY_ENERGY
	if(i - 1 >= 0) energy += dangle5_37[type][seq_int[i - 1]];
	if(j + 1 < seq_n) energy += dangle3_37[type][seq_int[j + 1]];
#else
	if(i - 1 >= 0 && j + 1 < seq_n) energy += mismatchExt37[type][seq_int[i - 1]][seq_int[j + 1]];
	else if(i - 1 >= 0) energy += dangle5_37[type][seq_int[i - 1]];
	else if(j + 1 < seq_n) energy += dangle3_37[type][seq_int[j + 1]];
#endif

	if(type > 2) energy += TerminalAU37;

	return energy;
}


// calc energy where bases in external [i, j] are unpaired
Float LinCapR::energy_external_unpaired(const int i, const int j) const{
	return 0;
}
