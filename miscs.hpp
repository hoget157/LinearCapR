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


const int BP_pair[NBASE][NBASE]=
/* @  A  C  G  U*/
{{ 0, 0, 0, 0, 0},
 { 0, 0, 0, 0, 5},
 { 0, 0, 0, 1, 0},
 { 0, 0, 2, 0, 3},
 { 0, 6, 0, 4, 0}};


/** Infinity as used in minimization routines */
#define INF 10000000


// convert base character to number
inline int base_to_num(const char base){
	if(base == 'A' || base == 'a') return 1;
	if(base == 'C' || base == 'c') return 2;
	if(base == 'G' || base == 'g') return 3;
	if(base == 'T' || base == 't' || base == 'U' || base == 'u') return 4;
	return 0;
}


// log space: borrowed from CONTRAfold
inline Float Fast_LogExpPlusOne(Float x){
  
    // Bounds for tolerance of 7.05e-06: (0, 11.8625)
    // Approximating interval: (0, 0.661537) --> ((T(-0.0065591595)*x+T(0.1276442762))*x+T(0.4996554598))*x+T(0.6931542306);
    // Approximating interval: (0.661537, 1.63202) --> ((T(-0.0155157557)*x+T(0.1446775699))*x+T(0.4882939746))*x+T(0.6958092989);
    // Approximating interval: (1.63202, 2.49126) --> ((T(-0.0128909247)*x+T(0.1301028251))*x+T(0.5150398748))*x+T(0.6795585882);
    // Approximating interval: (2.49126, 3.37925) --> ((T(-0.0072142647)*x+T(0.0877540853))*x+T(0.6208708362))*x+T(0.5909675829);
    // Approximating interval: (3.37925, 4.42617) --> ((T(-0.0031455354)*x+T(0.0467229449))*x+T(0.7592532310))*x+T(0.4348794399);
    // Approximating interval: (4.42617, 5.78907) --> ((T(-0.0010110698)*x+T(0.0185943421))*x+T(0.8831730747))*x+T(0.2523695427);
    // Approximating interval: (5.78907, 7.81627) --> ((T(-0.0001962780)*x+T(0.0046084408))*x+T(0.9634431978))*x+T(0.0983148903);
    // Approximating interval: (7.81627, 11.8625) --> ((T(-0.0000113994)*x+T(0.0003734731))*x+T(0.9959107193))*x+T(0.0149855051);
    // 8 polynomials needed.
    
    // assert(Float(0.0000000000) <= x && x <= Float(11.8624794162) && "Argument out-of-range.");
    if (x < Float(3.3792499610))
    {
        if (x < Float(1.6320158198))
        {
            if (x < Float(0.6615367791))
                return ((Float(-0.0065591595)*x+Float(0.1276442762))*x+Float(0.4996554598))*x+Float(0.6931542306);
            return ((Float(-0.0155157557)*x+Float(0.1446775699))*x+Float(0.4882939746))*x+Float(0.6958092989);
        }
        if (x < Float(2.4912588184))
            return ((Float(-0.0128909247)*x+Float(0.1301028251))*x+Float(0.5150398748))*x+Float(0.6795585882);
        return ((Float(-0.0072142647)*x+Float(0.0877540853))*x+Float(0.6208708362))*x+Float(0.5909675829);
    }
    if (x < Float(5.7890710412))
    {
        if (x < Float(4.4261691294))
            return ((Float(-0.0031455354)*x+Float(0.0467229449))*x+Float(0.7592532310))*x+Float(0.4348794399);
        return ((Float(-0.0010110698)*x+Float(0.0185943421))*x+Float(0.8831730747))*x+Float(0.2523695427);
    }
    if (x < Float(7.8162726752))
        return ((Float(-0.0001962780)*x+Float(0.0046084408))*x+Float(0.9634431978))*x+Float(0.0983148903);
    return ((Float(-0.0000113994)*x+Float(0.0003734731))*x+Float(0.9959107193))*x+Float(0.0149855051);
}


// returns z; e^z = e^x + e^y
inline Float logsumexp(Float x, Float y){
	if(x < y) swap(x, y);
	return (x - y < Float(11.8624794162) ? y + Fast_LogExpPlusOne(x - y) : x);
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
}


// apply effects of add_range()
inline void prefix_sum(vector<Float> &v){
	for(int i = 1; i < (int)v.size(); i++) v[i] += v[i - 1];
}


// returns whether [i, j] in t
inline bool contains(const Table &t, const int i, const int j){
	return t[j].count(i);
}
