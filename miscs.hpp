#pragma once

#include "utils.hpp"

#include <cmath>

/** The number of distinguishable base pairs */
#define NBPAIRS 7


const int BP_pair[5][5]=
/* @  A  C  G  U*/
{{ 0, 0, 0, 0, 0},
 { 0, 0, 0, 0, 5},
 { 0, 0, 0, 1, 0},
 { 0, 0, 2, 0, 3},
 { 0, 6, 0, 4, 0}};


/* rtype[pair[i][j]]:=pair[j][i] */
const int rtype[7] = {0, 2, 1, 4, 3, 6, 5};


/** Infinity as used in minimization routines */
#define INF 10000000


// 
inline int base_to_num(const char base){
	if(base == 'A' || base == 'a') return 1;
	if(base == 'C' || base == 'c') return 2;
	if(base == 'G' || base == 'g') return 3;
	if(base == 'T' || base == 't' || base == 'U' || base == 'u') return 4;
	return 0;
}


// returns z; e^z = e^x + e^y
inline Float logsumexp(const Float x, const Float y){
	if(x == -INF) return y;
	if(y == -INF) return x;
	return (x > y ? x + log1p(exp(y - x)) : y + log1p(exp(x - y)));
}


// t[i, j] += score
inline Float update_sum(Table &t, const int i, const int j, const Float score){
	// DUMP("update", i, j, score);
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


	// inline Float logsumexp_equal(Float &x, const Float y) const{
	// 	return x = logsumexp(x, y);
	// }

