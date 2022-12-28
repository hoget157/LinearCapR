#pragma once

#include "types.hpp"

#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cmath>


// dump vectors
template<class T>
ostream &operator<<(ostream &os,const vector<T> &v){
	// os << fixed;
	for(int i = 0;i < (int)v.size();i++) os << (i ? "," : "[") /*<< setprecision(15)*/ << v[i];
	os << "]";
	return os;
}


// dump unordered_map
template<class T1, class T2>
ostream &operator<<(ostream &os,const unordered_map<T1, T2> &m){
	os << "[";
	for(const auto [x, y] : m) os << "(" << x << ", " << y /** -kT*/ << "), ";
	os << "]";
	return os;
}


// update min/max
template<class T> inline T &chmin(T &a,const T &b){ return a = min(a,b); }
template<class T> inline T &chmax(T &a,const T &b){ return a = max(a,b); }


// dump variables
void DUMP();
template <class THead, class... TTail>
void DUMP(THead&& head, TTail&&... tail){
	cout << head << ", ";
	DUMP(move(tail)...);
}


// class to measure time
template<class T>
struct Time{
	chrono::system_clock::time_point begin_time;

	void init(){
		begin_time = chrono::system_clock::now();
	}

	double measure() const{
		return (double)chrono::duration_cast<T>(chrono::system_clock::now() - begin_time).count();
	}
};
