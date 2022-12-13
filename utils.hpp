#pragma once

#include <iostream>
#include <vector>
#include <chrono>

using namespace std;

// dump vectors
template<class T>
ostream &operator<<(ostream &os,const vector<T> &v){
	for(int i = 0;i < (int)v.size();i++) os << (i ? "," : "[") << v[i];
	os << "]";
	return os;
}

// update min/max
template<class T> T &chmin(T &a,const T &b){ return a = min(a,b); }
template<class T> T &chmax(T &a,const T &b){ return a = max(a,b); }

// class to measure time
template<class T>
struct Time{
	chrono::system_clock::time_point begin_time;

	void init(){
		begin_time = chrono::system_clock::now();
	}

	double measure(){
		return (double)chrono::duration_cast<T>(chrono::system_clock::now() - begin_time).count();
	}
};
