#include "LinCapR.hpp"
#include "beam_prune.hpp"
#include "dp_table_api.hpp"
#include "seq_utils.hpp"
#include "energy_linearcapr.hpp"
#include "energy_raccess.hpp"

#include <fstream>
#include <algorithm>
#include <iostream>
#include <cmath>

LinCapR::LinCapR(int beam_size, energy::Model model, EnergyEngine engine, bool normalize_profiles, Float normalize_warn_eps)
	: params(energy::get_params(model)),
	  beam_size(beam_size),
	  normalize_profiles(normalize_profiles),
	  normalize_warn_eps(normalize_warn_eps) {
	switch (engine) {
	case EnergyEngine::Raccess:
		_energy = lcr::make_raccess_energy();
		break;
	case EnergyEngine::LinearCapR:
	default:
		_energy.reset(new lcr::LinearCapREnergyModel(params));
		break;
	}
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
	const Float logZ = alpha_O[seq_n - 1];
	if (dynamic_cast<const lcr::LinearCapREnergyModel*>(_energy.get()) != nullptr) {
		return (logZ * -(params.temperature + params.k0) * params.gas_constant) / 1000;
	}
	return -_energy->kT() * logZ;
}

Float LinCapR::get_logZ() const{
	return alpha_O[seq_n - 1];
}

const vector<Float>& LinCapR::get_prob_stem() const{
	return prob_S;
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
	_energy->set_sequence(seq, seq_int);

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
				lcr::dp::update_sum(alpha_S, i - 1, j + 1, score - _energy->energy_loop(i - 1, j + 1, i, j) / _energy->kT());
			}
			
			// M2 -> S
			for(int n = 0; n <= MULTI_MAX_UNPAIRED && j + n < seq_n; n++){
				lcr::dp::update_sum(alpha_M2, i, j + n, score - (_energy->energy_multi_bif(i, j) + _energy->energy_multi_unpaired(j + 1, j + n)) / _energy->kT());
			}

			// SE -> S: p..i..j..q, [p - 1, q] can be pair
			for(int p = i; i - p <= MAXLOOP && p >= 1; p--){
				for(int q = next_pair[seq_int[p - 1]][j + 1]; q < seq_n && (q - j - 1) + (i - p) <= MAXLOOP; q = next_pair[seq_int[p - 1]][q + 1]){
					if((p == i && q == j + 1)) continue;
					lcr::dp::update_sum(alpha_SE, p, q - 1, score - _energy->energy_loop(p - 1, q, i, j) / _energy->kT());
				}
			}

			// O -> O + S
			lcr::dp::update_sum(alpha_O, j, (i - 1 >= 0 ? alpha_O[i - 1] : 0) + score - _energy->energy_external(i, j) / _energy->kT());
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
				lcr::dp::update_sum(alpha_SE, i, j, score - _energy->energy_multi_closing(i - 1, j + 1) / _energy->kT());
			}
		}

		// SE -> (Hairpin)
		for(int n = TURN; n <= MAXLOOP; n++){
			const int i = j - n + 1;
			if(i - 1 >= 0 && j + 1 < seq_n && can_pair(i - 1, j + 1)){
				lcr::dp::update_sum(alpha_SE, i, j, -_energy->energy_hairpin(i - 1, j + 1) / _energy->kT());
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
			lcr::dp::update_sum(alpha_O, j + 1, alpha_O[j] - _energy->energy_external_unpaired(j + 1, j + 1) / _energy->kT());
		}
	}
}


// calc outside variables
void LinCapR::calc_outside(){
	for(int j = seq_n - 1; j >= 0; j--){
		// O
		// O -> O
		lcr::dp::update_sum(beta_O, j, (j + 1 < seq_n ? beta_O[j + 1] : 0) - _energy->energy_external_unpaired(j + 1, j + 1) / _energy->kT());
		
		// O -> O + S
		for(const auto [i, score] : alpha_S[j]){
			lcr::dp::update_sum(beta_O, i, score + (j + 1 < seq_n ? beta_O[j + 1] : 0) - _energy->energy_external(i, j) / _energy->kT());
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
				lcr::dp::update_sum(beta_M, i, j, lcr::dp::get_value(beta_SE, i, j) - _energy->energy_multi_closing(i - 1, j + 1) / _energy->kT());
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
			lcr::dp::update_sum(beta_S, i, j, (i - 1 >= 0 ? alpha_O[i - 1] : 0) + (j + 1 < seq_n ? beta_O[j + 1] : 0) - _energy->energy_external(i, j) / _energy->kT());

			// SE -> S
			for(int p = i; i - p <= MAXLOOP && p >= 1; p--){
				for(int q = next_pair[seq_int[p - 1]][j + 1]; q < seq_n && (q - j - 1) + (i - p) <= MAXLOOP; q = next_pair[seq_int[p - 1]][q + 1]){
					if((p == i && q == j + 1)) continue;
					lcr::dp::update_sum(beta_S, i, j, lcr::dp::get_value(beta_SE, p, q - 1) - _energy->energy_loop(p - 1, q, i, j) / _energy->kT());
				}
			}

			// S -> S
			if(i - 1 >= 0 && j + 1 < seq_n){
				lcr::dp::update_sum(beta_S, i, j, lcr::dp::get_value(beta_S, i - 1, j + 1) - _energy->energy_loop(i - 1, j + 1, i, j) / _energy->kT());
			}
			
			// M2 -> S
			for(int n = 0; n <= MULTI_MAX_UNPAIRED && j + n < seq_n; n++){
				lcr::dp::update_sum(beta_S, i, j, lcr::dp::get_value(beta_M2, i, j + n) - (_energy->energy_multi_bif(i, j) + _energy->energy_multi_unpaired(j + 1, j + n)) / _energy->kT());
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
			lcr::dp::add_range(prob_H, j, k, exp(score - _energy->energy_hairpin(j - 1, k + 1) / _energy->kT() - logZ));

			// B, I
			for(int p = j; p <= min(j + MAXLOOP, k - 1); p++){
				for(int q = k; q >= p + TURN + 1 && (p - j) + (k - q) <= MAXLOOP; q--){
					if((p == j && q == k) || !lcr::dp::contains(alpha_S, p, q)) continue;
					auto it_a = alpha_S[q].find(p);
					if(it_a == alpha_S[q].end()) continue;
					const Float new_score = exp(score + it_a->second - _energy->energy_loop(j - 1, k + 1, p, q) / _energy->kT() - logZ);
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
				const Float new_score = exp(score + beta_M[k][j] - _energy->energy_multi_unpaired(j, p - 1) / _energy->kT() - logZ);
				lcr::dp::add_range(prob_M, j, p - 1, new_score);
			}
		}
	}
	for(int q = 0; q < seq_n; q++){
		for(const auto [j, score] : alpha_S[q]){
			for(int k = q + 1; k <= min(seq_n - 1, q + MAXLOOP); k++){
				if(!lcr::dp::contains(beta_M2, j, k)) continue;
				const Float new_score = exp(score + beta_M2[k][j] - (_energy->energy_multi_bif(j, q) + _energy->energy_multi_unpaired(q + 1, k)) / _energy->kT() - logZ);
				lcr::dp::add_range(prob_M, q + 1, k, new_score);
			}
		}
	}
	lcr::dp::prefix_sum(prob_M);

	// S
	for(int j = 0; j < seq_n; j++){
		for(const auto [i, score] : alpha_S[j]){
			auto it_b = beta_S[j].find(i);
			if(it_b == beta_S[j].end()) continue;
			const Float new_score = exp(score + it_b->second - logZ);
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
		if(fabs(sum_prob_i - 1.0) > normalize_warn_eps){
			cerr << "warn: prob_sum[" << i << "]=" << sum_prob_i << endl;
		}
		if(!normalize_profiles) continue;
		// sum of probabilities to 1
		for(int j = 0; j < NPROBS; j++) probs[j]->at(i) /= sum_prob_i;
	}
}

void LinCapR::debug_stem_pairs(int idx, int topn) const {
	if(idx < 0 || idx >= seq_n){
		cerr << "debug_stem_pairs: idx out of range: " << idx << endl;
		return;
	}
	if(topn <= 0) topn = 1;

	struct Item {
		int i;
		int j;
		double prob;
	};
	vector<Item> items;
	items.reserve(seq_n);

	const double logZ = alpha_O[seq_n - 1];
	double total = 0.0;
	int missing_beta = 0;

	// idx as left endpoint (i = idx, j varies)
	for(int j = 0; j < seq_n; j++){
		auto it_a = alpha_S[j].find(idx);
		if(it_a == alpha_S[j].end()) continue;
		auto it_b = beta_S[j].find(idx);
		if(it_b == beta_S[j].end()){
			missing_beta++;
			continue;
		}
		const double prob = exp(it_a->second + it_b->second - logZ);
		if(prob <= 0.0) continue;
		items.push_back({idx, j, prob});
		total += prob;
	}

	// idx as right endpoint (j = idx, i varies)
	for(const auto& kv : alpha_S[idx]){
		const int i = kv.first;
		auto it_b = beta_S[idx].find(i);
		if(it_b == beta_S[idx].end()){
			missing_beta++;
			continue;
		}
		const double prob = exp(kv.second + it_b->second - logZ);
		if(prob <= 0.0) continue;
		items.push_back({i, idx, prob});
		total += prob;
	}

	sort(items.begin(), items.end(), [](const Item& a, const Item& b){
		return a.prob > b.prob;
	});

	cerr << "debug_stem_pairs i=" << idx
	     << " total=" << total
	     << " prob_S=" << prob_S[idx]
	     << " missing_beta=" << missing_beta << endl;
	const int limit = min(topn, static_cast<int>(items.size()));
	for(int k = 0; k < limit; k++){
		cerr << "  pair (" << items[k].i << "," << items[k].j << ") prob=" << items[k].prob << endl;
	}
}

void LinCapR::debug_pair(int i, int j) const {
	if(i < 0 || j < 0 || i >= seq_n || j >= seq_n || i >= j){
		cerr << "debug_pair: invalid i,j: " << i << "," << j << endl;
		return;
	}
	const double logZ = alpha_O[seq_n - 1];
	const bool pairable = can_pair(i, j);
	cerr << "debug_pair (" << i << "," << j << ")"
	     << " bases=" << seq[i] << "," << seq[j]
	     << " can_pair=" << pairable << endl;

	auto it_a = alpha_S[j].find(i);
	auto it_b = beta_S[j].find(i);
	if(it_a == alpha_S[j].end()){
		cerr << "  alpha_S missing" << endl;
	} else {
		cerr << "  alpha_S=" << it_a->second << endl;
	}
	if(it_b == beta_S[j].end()){
		cerr << "  beta_S missing" << endl;
	} else {
		cerr << "  beta_S=" << it_b->second << endl;
	}
	if(it_a != alpha_S[j].end() && it_b != beta_S[j].end()){
		const double prob = exp(it_a->second + it_b->second - logZ);
		cerr << "  pair_prob=" << prob << endl;
	}

	if(i - 1 >= 0 && j + 1 < seq_n){
		const double loop_energy = _energy->energy_loop(i - 1, j + 1, i, j);
		cerr << "  energy_loop(i-1,j+1,i,j)=" << loop_energy << endl;
	}
}

void LinCapR::debug_prob(int idx) const {
	if(idx < 0 || idx >= seq_n){
		cerr << "debug_prob: idx out of range: " << idx << endl;
		return;
	}
	const double b = prob_B[idx];
	const double h = prob_H[idx];
	const double in = prob_I[idx];
	const double m = prob_M[idx];
	const double e = prob_E[idx];
	const double s = prob_S[idx];
	const double sum = b + h + in + m + e + s;
	cerr << "debug_prob i=" << idx
	     << " B=" << b
	     << " H=" << h
	     << " I=" << in
	     << " M=" << m
	     << " E=" << e
	     << " S=" << s
	     << " sum=" << sum << endl;
}

void LinCapR::debug_internal(int idx, int topn) const {
	if(idx < 0 || idx >= seq_n){
		cerr << "debug_internal: idx out of range: " << idx << endl;
		return;
	}
	if(topn <= 0) topn = 1;

	struct Item {
		int j;
		int k;
		int p;
		int q;
		double prob;
		const char* side;
	};
	vector<Item> items;
	const double logZ = alpha_O[seq_n - 1];

	for(int k = 0; k < seq_n; k++){
		for(const auto [j, score] : beta_SE[k]){
			for(int p = j; p <= min(j + MAXLOOP, k - 1); p++){
				for(int q = k; q >= p + TURN + 1 && (p - j) + (k - q) <= MAXLOOP; q--){
					if((p == j && q == k) || !lcr::dp::contains(alpha_S, p, q)) continue;
					auto it_a = alpha_S[q].find(p);
					if(it_a == alpha_S[q].end()) continue;
					const Float new_score = exp(score + it_a->second - _energy->energy_loop(j - 1, k + 1, p, q) / _energy->kT() - logZ);

					const bool left_internal = (q != k) && (idx >= j && idx <= (p - 1));
					const bool right_internal = (p != j) && (idx >= (q + 1) && idx <= k);
					if(left_internal){
						items.push_back({j, k, p, q, new_score, "left"});
					}
					if(right_internal){
						items.push_back({j, k, p, q, new_score, "right"});
					}
				}
			}
		}
	}

	sort(items.begin(), items.end(), [](const Item& a, const Item& b){
		return a.prob > b.prob;
	});

	const auto base_at = [&](int pos) -> char {
		if(pos < 0 || pos >= seq_n) return 'N';
		return seq[pos];
	};

	cerr << "debug_internal i=" << idx << " candidates=" << items.size() << endl;
	const int limit = min(topn, static_cast<int>(items.size()));
	for(int t = 0; t < limit; t++){
		const auto& it = items[t];
		const int outer_i = it.j - 1;
		const int outer_j = it.k + 1;
		const int left_len = it.p - it.j;
		const int right_len = it.k - it.q;
		const bool outer_pairable = (0 <= outer_i && outer_j < seq_n) ? can_pair(outer_i, outer_j) : false;
		const bool outer_in_alpha = (0 <= outer_i && outer_j < seq_n) ? lcr::dp::contains(alpha_S, outer_i, outer_j) : false;
		const double loop_energy = (outer_i >= 0 && outer_j < seq_n)
			? _energy->energy_loop(outer_i, outer_j, it.p, it.q)
			: 0.0;
		cerr << "  " << it.side
		     << " (j,k,p,q)=(" << it.j << "," << it.k << "," << it.p << "," << it.q << ")"
		     << " prob=" << it.prob
		     << " outer=(" << outer_i << "," << outer_j << ")"
		     << " inner=(" << it.p << "," << it.q << ")"
		     << " len=(" << left_len << "," << right_len << ")"
		     << " bases outer=" << base_at(outer_i) << "," << base_at(outer_j)
		     << " inner=" << base_at(it.p) << "," << base_at(it.q)
		     << " outer_can_pair=" << outer_pairable
		     << " outer_in_alpha=" << outer_in_alpha
		     << " loopE=" << loop_energy
		     << endl;
	}
}

void LinCapR::debug_multi_unpaired(int i, int j) const {
	if(i < 0 || j < 0 || i >= seq_n || j >= seq_n || i > j){
		cerr << "debug_multi_unpaired: invalid range: " << i << "," << j << endl;
		return;
	}
	const int len_closed = j - i + 1;
	const int len_half = j - i;
	const double energy = _energy->energy_multi_unpaired(i, j);
	cerr << "debug_multi_unpaired (i,j)=(" << i << "," << j << ")"
	     << " len_closed=" << len_closed
	     << " len_half=" << len_half
	     << " energy=" << energy
	     << endl;
}

void LinCapR::debug_multi_prob(int idx, int topn) const {
	if(idx < 0 || idx >= seq_n){
		cerr << "debug_multi_prob: idx out of range: " << idx << endl;
		return;
	}
	if(topn <= 0) topn = 1;
	const Float logZ = alpha_O[seq_n - 1];
	double from_mb = 0.0;
	double from_s = 0.0;

	struct Item {
		int j;
		int q;
		int k;
		double prob;
	};
	vector<Item> items;

	// Contribution from alpha_MB / beta_M (calc_profile M part 1)
	for(int k = 0; k < seq_n; k++){
		for(const auto [p, score] : alpha_MB[k]){
			const int j_min = max(0, p - MAXLOOP);
			for(int j = p - 1; j >= j_min; j--){
				if(!lcr::dp::contains(beta_M, j, k)) continue;
				if(idx < j || idx > (p - 1)) continue;
				auto it_b = beta_M[k].find(j);
				if(it_b == beta_M[k].end()) continue;
				const Float new_score = exp(score + it_b->second
					- _energy->energy_multi_unpaired(j, p - 1) / _energy->kT()
					- logZ);
				from_mb += new_score;
			}
		}
	}

	// Contribution from alpha_S / beta_M2 (calc_profile M part 2)
	for(int q = 0; q < seq_n; q++){
		for(const auto [j, score] : alpha_S[q]){
			const int k_max = min(seq_n - 1, q + MAXLOOP);
			for(int k = q + 1; k <= k_max; k++){
				if(!lcr::dp::contains(beta_M2, j, k)) continue;
				if(idx < (q + 1) || idx > k) continue;
				auto it_b = beta_M2[k].find(j);
				if(it_b == beta_M2[k].end()) continue;
				const Float new_score = exp(score + it_b->second
					- (_energy->energy_multi_bif(j, q) + _energy->energy_multi_unpaired(q + 1, k)) / _energy->kT()
					- logZ);
				from_s += new_score;
				items.push_back({j, q, k, new_score});
			}
		}
	}

	const double total = from_mb + from_s;
	sort(items.begin(), items.end(), [](const Item& a, const Item& b){
		return a.prob > b.prob;
	});

	cerr << "debug_multi_prob i=" << idx
	     << " M_total=" << prob_M[idx]
	     << " M_from_mb=" << from_mb
	     << " M_from_s=" << from_s
	     << " M_sum=" << total
	     << endl;

	const int limit = min(topn, static_cast<int>(items.size()));
	for(int t = 0; t < limit; t++){
		const auto& it = items[t];
		cerr << "  M2 (j,q,k)=(" << it.j << "," << it.q << "," << it.k << ")"
		     << " prob=" << it.prob
		     << " range=(" << (it.q + 1) << "," << it.k << ")"
		     << endl;
	}
}


// calc energy of hairpin loop [i, j]
