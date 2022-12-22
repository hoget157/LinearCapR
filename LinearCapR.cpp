#include "LinearCapR.hpp"
#include "utils.hpp"

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
Float LinearCapR::prune(unordered_map<int, Float> &states){
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
}
