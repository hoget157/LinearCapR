#pragma once

#include <cmath>
#include <vector>
#include <unordered_map>

using namespace std;


// using Float = float;
using Float = double;
// using Float = long double;

using Table = vector<unordered_map<int, Float>>;


/** The number of distinguishable base pairs */
#define NBPAIRS 7
/** The number of distinguishable base*/
#define NBASE 5
/** The number of table*/
#define NTABLES 6
/** The number of type of profiles*/
#define NPROBS 6
/** Infinity as used in minimization routines */
#define INF 10000000


// convert bases to base pair index
const int BP_pair[NBASE][NBASE]=
/* @  A  C  G  U*/
{{ 0, 0, 0, 0, 0},
 { 0, 0, 0, 0, 5},
 { 0, 0, 0, 1, 0},
 { 0, 0, 2, 0, 3},
 { 0, 6, 0, 4, 0}};


// convert base character to number
inline int base_to_num(const char base){
	if(base == 'A' || base == 'a') return 1;
	if(base == 'C' || base == 'c') return 2;
	if(base == 'G' || base == 'g') return 3;
	if(base == 'T' || base == 't' || base == 'U' || base == 'u') return 4;
	return 0;
}


inline Float logsumexp(const Float x, const Float y){
	if(x == -INF) return y;
	if(y == -INF) return x;
	return (x > y ? x + log1p(exp(y - x)) : y + log1p(exp(x - y)));
}


// t[i, j] += score
inline Float update_sum(Table &t, const int i, const int j, const Float score){
	return t[j][i] = (t[j].count(i) ? logsumexp(t[j][i], score) : score);
}


// v[i] += score
inline Float update_sum(vector<Float> &v, const int i, const Float score){
	return v[i] = logsumexp(v[i], score);
}


// returns t[i, j] if exists, else default value
inline Float get_value(const Table &t, const int i, const int j, const Float default_value = -INF){
	return (t[j].count(i) ? t[j].at(i) : default_value);
}


// for k in [i, j]: v[k] += x
// call prefix_sum() to complete
inline void add_range(vector<Float> &v, const int i, const int j, const Float x){
	v[i] += x;
	if(j + 1 < (int)v.size()) v[j + 1] -= x;
    // for(int k = i; k <= j; k++) v[k] += x;
}


// apply effects of add_range()
inline void prefix_sum(vector<Float> &v){
	for(int i = 1; i < (int)v.size(); i++) v[i] += v[i - 1];
}


// returns whether [i, j] in t
inline bool contains(const Table &t, const int i, const int j){
	return t[j].count(i);
}
