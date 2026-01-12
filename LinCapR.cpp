#include "LinCapR.hpp"
#include "beam_prune.hpp"
#include "dp_table_api.hpp"
#include "seq_utils.hpp"

#include <fstream>
#include <algorithm>
#include <cstring>

LinCapR::LinCapR(int beam_size, energy::Model model)
	: params(energy::get_params(model)), beam_size(beam_size){
	if(params.use_fast_logsumexp) set_logsumexp_fast_mode();
	else set_logsumexp_legacy_mode();
}


// prune top-k states
Float LinCapR::prune(Map<int, Float> &states) const{
	return lcr::beam::prune_states(states, beam_size,
				       [this](const int i, const Float score) {
					 return (i >= 1 ? alpha_O[i - 1] : Float(0)) + score;
				       });
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

	ofs << "Multiloop ";
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
		seq_int[i] = lcr::seq::base_to_num(seq[i]);
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
	lcr::seq::build_next_pair(seq_int, next_pair);
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
				lcr::dp::update_sum(alpha_S, i - 1, j + 1, score - energy_loop(i - 1, j + 1, i, j) / params.kT);
			}
			
			// M2 -> S
			for(int n = 0; n <= MULTI_MAX_UNPAIRED && j + n < seq_n; n++){
				lcr::dp::update_sum(alpha_M2, i, j + n, score - (energy_multi_bif(i, j) + energy_multi_unpaired(j + 1, j + n)) / params.kT);
			}

			// SE -> S: p..i..j..q, [p - 1, q] can be pair
			for(int p = i; i - p <= MAXLOOP && p >= 1; p--){
				for(int q = next_pair[seq_int[p - 1]][j + 1]; q < seq_n && (q - j - 1) + (i - p) <= MAXLOOP; q = next_pair[seq_int[p - 1]][q + 1]){
					if((p == i && q == j + 1)) continue;
					lcr::dp::update_sum(alpha_SE, p, q - 1, score - energy_loop(p - 1, q, i, j) / params.kT);
				}
			}

			// O -> O + S
			lcr::dp::update_sum(alpha_O, j, (i - 1 >= 0 ? alpha_O[i - 1] : 0) + score - energy_external(i, j) / params.kT);
		}

		// M2
		prune(alpha_M2[j]);
		for(const auto [i, score] : alpha_M2[j]){
			// M1 -> M2
			lcr::dp::update_sum(alpha_M1, i, j, score);

			// MB -> M1 + M2
			if(i - 1 >= 0){
				for(const auto [k, score_m1] : alpha_M1[i - 1]){
					lcr::dp::update_sum(alpha_MB, k, j, score_m1 + score);
				}
			}
		}

		// MB
		prune(alpha_MB[j]);
		for(const auto [i, score] : alpha_MB[j]){
			// M1 -> MB
			lcr::dp::update_sum(alpha_M1, i, j, score);

			// M -> MB
			for(int n = 0; n <= MULTI_MAX_UNPAIRED && i - n >= 0; n++){
				lcr::dp::update_sum(alpha_M, i - n, j, score);
			}
		}

		// M1
		prune(alpha_M1[j]);

		// M
		prune(alpha_M[j]);
		for(const auto [i, score] : alpha_M[j]){
			// SE -> M
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				lcr::dp::update_sum(alpha_SE, i, j, score - energy_multi_closing(i - 1, j + 1) / params.kT);
			}
		}

		// SE -> (Hairpin)
		for(int n = TURN; n <= MAXLOOP; n++){
			const int i = j - n + 1;
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				lcr::dp::update_sum(alpha_SE, i, j, -energy_hairpin(i - 1, j + 1) / params.kT);
			}
		}

		// SE
		prune(alpha_SE[j]);
		for(const auto [i, score] : alpha_SE[j]){
			// S -> SE
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				lcr::dp::update_sum(alpha_S, i - 1, j + 1, score);
			}
		}

		// O -> O
		if(j + 1 < seq_n){
			lcr::dp::update_sum(alpha_O, j + 1, alpha_O[j] - energy_external_unpaired(j + 1, j + 1) / params.kT);
		}
	}
}


// calc outside variables
void LinCapR::calc_outside(){
	for(int j = seq_n - 1; j >= 0; j--){
		// O
		// O -> O
		lcr::dp::update_sum(beta_O, j, (j + 1 < seq_n ? beta_O[j + 1] : 0) - energy_external_unpaired(j + 1, j + 1) / params.kT);
		
		// O -> O + S
		for(const auto [i, score] : alpha_S[j]){
			lcr::dp::update_sum(beta_O, i, score + (j + 1 < seq_n ? beta_O[j + 1] : 0) - energy_external(i, j) / params.kT);
		}

		// SE
		for(const auto [i, _] : alpha_SE[j]){
			// S -> SE
			if(i - 1 >= 0 && j + 1 < seq_n){
				lcr::dp::update_sum(beta_SE, i, j, lcr::dp::get_value(beta_S, i - 1, j + 1));
			}
		}

		// M
		for(const auto [i, _] : alpha_M[j]){
			// SE -> M
			if(i - 1 >= 0 && j + 1 < seq_n){
				lcr::dp::update_sum(beta_M, i, j, lcr::dp::get_value(beta_SE, i, j) - energy_multi_closing(i - 1, j + 1) / params.kT);
			}
		}

		// MB
		for(const auto [i, _] : alpha_MB[j]){
			// M1 -> MB
			lcr::dp::update_sum(beta_MB, i, j, lcr::dp::get_value(beta_M1, i, j));

			// M -> MB
			for(int n = 0; n <= MULTI_MAX_UNPAIRED && i - n >= 0; n++){
				lcr::dp::update_sum(beta_MB, i, j, lcr::dp::get_value(beta_M, i - n, j));
			}
		}

		// M1, M2
		for(const auto [i, score_M2] : alpha_M2[j]){
			// M1 -> M2
			lcr::dp::update_sum(beta_M2, i, j, lcr::dp::get_value(beta_M1, i, j));

			// MB -> M1 + M2
			if(i - 1 < 0) continue;
			for(const auto [k, score_M1] : alpha_M1[i - 1]){
				lcr::dp::update_sum(beta_M1, k, i - 1, lcr::dp::get_value(beta_MB, k, j) + score_M2);
				lcr::dp::update_sum(beta_M2, i, j, lcr::dp::get_value(beta_MB, k, j) + score_M1);
			}
		}

		// S
		for(const auto [i, _] : alpha_S[j]){
			// O -> O + S
			lcr::dp::update_sum(beta_S, i, j, (i - 1 >= 0 ? alpha_O[i - 1] : 0) + (j + 1 < seq_n ? beta_O[j + 1] : 0) - energy_external(i, j) / params.kT);

			// SE -> S
			for(int p = i; i - p <= MAXLOOP && p >= 1; p--){
				for(int q = next_pair[seq_int[p - 1]][j + 1]; q < seq_n && (q - j - 1) + (i - p) <= MAXLOOP; q = next_pair[seq_int[p - 1]][q + 1]){
					if((p == i && q == j + 1)) continue;
					lcr::dp::update_sum(beta_S, i, j, lcr::dp::get_value(beta_SE, p, q - 1) - energy_loop(p - 1, q, i, j) / params.kT);
				}
			}

			// S -> S
			if(i - 1 >= 0 && j + 1 < seq_n){
				lcr::dp::update_sum(beta_S, i, j, lcr::dp::get_value(beta_S, i - 1, j + 1) - energy_loop(i - 1, j + 1, i, j) / params.kT);
			}
			
			// M2 -> S
			for(int n = 0; n <= MULTI_MAX_UNPAIRED && j + n < seq_n; n++){
				lcr::dp::update_sum(beta_S, i, j, lcr::dp::get_value(beta_M2, i, j + n) - (energy_multi_bif(i, j) + energy_multi_unpaired(j + 1, j + n)) / params.kT);
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
			lcr::dp::add_range(prob_H, j, k, exp(score - energy_hairpin(j - 1, k + 1) / params.kT - logZ));

			// B, I
			for(int p = j; p <= min(j + MAXLOOP, k - 1); p++){
				for(int q = k; q >= p + TURN + 1 && (p - j) + (k - q) <= MAXLOOP; q--){
					if((p == j && q == k) || !lcr::dp::contains(alpha_S, p, q)) continue;
					const Float new_score = exp(score + alpha_S[q][p] - energy_loop(j - 1, k + 1, p, q) / params.kT - logZ);
					lcr::dp::add_range((q == k ? prob_B : prob_I), j, p - 1, new_score);
					lcr::dp::add_range((p == j ? prob_B : prob_I), q + 1, k, new_score);
				}
			}
		}
	}
	lcr::dp::prefix_sum(prob_B);
	lcr::dp::prefix_sum(prob_H);
	lcr::dp::prefix_sum(prob_I);

	// M
	for(int k = 0; k < seq_n; k++){
		for(const auto [p, score] : alpha_MB[k]){
			for(int j = p - 1; j >= max(0, p - MAXLOOP); j--){
				if(!lcr::dp::contains(beta_M, j, k)) continue;
				const Float new_score = exp(score + beta_M[k][j] - energy_multi_unpaired(j, p - 1) / params.kT - logZ);
				lcr::dp::add_range(prob_M, j, p - 1, new_score);
			}
		}
	}
	for(int q = 0; q < seq_n; q++){
		for(const auto [j, score] : alpha_S[q]){
			for(int k = q + 1; k <= min(seq_n - 1, q + MAXLOOP); k++){
				if(!lcr::dp::contains(beta_M2, j, k)) continue;
				const Float new_score = exp(score + beta_M2[k][j] - (energy_multi_bif(j, q) + energy_multi_unpaired(q + 1, k)) / params.kT - logZ);
				lcr::dp::add_range(prob_M, q + 1, k, new_score);
			}
		}
	}
	lcr::dp::prefix_sum(prob_M);

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
	if(!params.has_special_hairpins) return -1;

	const int d = j - i - 1;
	const char *loops_seq = nullptr;
	if(d == 3) loops_seq = params.Triloops;
	else if(d == 4) loops_seq = params.Tetraloops;
	else if(d == 6) loops_seq = params.Hexaloops;
	else return -1;

	if(!loops_seq) return -1;

	const string loop = seq.substr(i, d + 2);
	const char *sp = strstr(loops_seq, loop.c_str());
	return (sp ? (sp - loops_seq) / (d + 3) : -1);
}


// calc energy of hairpin loop [i, j]
Float LinCapR::energy_hairpin(const int i, const int j) const{
	const int type = BP_pair[seq_int[i]][seq_int[j]];
	const int d = j - i - 1;
	
	// check special hairpin
	const int index = special_hairpin(i, j);
	if(index != -1){
		if(d == 3 && params.Triloop37) return params.Triloop37[index];
		if(d == 4 && params.Tetraloop37) return params.Tetraloop37[index];
		if(d == 6 && params.Hexaloop37) return params.Hexaloop37[index];
	}

	// initiation
	Float energy = (d <= MAXLOOP ? params.hairpin37[d] : params.hairpin37[30] + params.lxc37 * log(d / 30.));
	
	if(d != 3){
		energy += params.mismatchH37[type][seq_int[i + 1]][seq_int[j - 1]];
	}else if(type > 2){
		energy += params.TerminalAU37;
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
		return params.stack37[type1][type2];
	}

	if(dmin == 0){
		// bulge
		Float energy = (d <= MAXLOOP ? params.bulge37[d] : params.bulge37[30] + params.lxc37 * log(d / 30.));

		if(dmax == 1) energy += params.stack37[type1][type2];
		else{
			if(type1 > 2) energy += params.TerminalAU37;
			if(type2 > 2) energy += params.TerminalAU37;
		}
		return energy;
	}

	// internal
	// specieal internal loops
	if(d1 == 1 && d2 == 1) return (*params.int11_37)[type1][type2][si][sj];
	if(d1 == 1 && d2 == 2) return (*params.int21_37)[type1][type2][si][sq][sj];
	if(d1 == 2 && d2 == 1) return (*params.int21_37)[type2][type1][sq][si][sp];
	if(d1 == 2 && d2 == 2) return (*params.int22_37)[type1][type2][si][sp][sq][sj];

	// generic internal loop
	Float energy =  (d <= MAXLOOP ? params.internal_loop37[d] : params.internal_loop37[30] + params.lxc37 * log(d / 30.));
	energy += min(params.MAX_NINIO, params.ninio37 * (dmax - dmin));
	
	// mismatch: different for sizes
	if(dmin == 1){ // 1xn
		energy += params.mismatch1nI37[type1][si][sj] + params.mismatch1nI37[type2][sq][sp];
	}else if(dmin == 2 && dmax == 3){ // 2x3
		energy += params.mismatch23I37[type1][si][sj] + params.mismatch23I37[type2][sq][sp];
	}else{ // others
		energy += params.mismatchI37[type1][si][sj] + params.mismatchI37[type2][sq][sp];
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
	return energy_multi_bif(j, i) + params.ML_closing37;
}


// calc energy of bifurcation [i, j] in a multiloop
Float LinCapR::energy_multi_bif(const int i, const int j) const{
	const int type = BP_pair[seq_int[i]][seq_int[j]];
	Float energy = params.ML_intern37;

	const bool has_left = (i - 1) >= 0;
	const bool has_right = (j + 1) < seq_n;

	if(params.allow_mismatch_multi && params.mismatchM37 && has_left && has_right){
		energy += params.mismatchM37[type][seq_int[i - 1]][seq_int[j + 1]];
	}else{
		if(has_left) energy += params.dangle5_37[type][seq_int[i - 1]];
		if(has_right) energy += params.dangle3_37[type][seq_int[j + 1]];
	}

	if(type > 2) energy += params.TerminalAU37;

	return energy;
}


// calc energy of external loop [i, j]
Float LinCapR::energy_external(const int i, const int j) const{
	const int type = BP_pair[seq_int[i]][seq_int[j]];
	Float energy = 0;

	const bool has_left = (i - 1) >= 0;
	const bool has_right = (j + 1) < seq_n;

	if(params.allow_mismatch_external && params.mismatchExt37 && has_left && has_right){
		energy += params.mismatchExt37[type][seq_int[i - 1]][seq_int[j + 1]];
	}else{
		if(has_left) energy += params.dangle5_37[type][seq_int[i - 1]];
		if(has_right) energy += params.dangle3_37[type][seq_int[j + 1]];
	}

	if(type > 2) energy += params.TerminalAU37;

	return energy;
}


// calc energy where bases in external [i, j] are unpaired
Float LinCapR::energy_external_unpaired(const int i, const int j) const{
	return 0;
}
