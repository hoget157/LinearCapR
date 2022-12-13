#pragma once

#include <string>
#include <vector>

using namespace std;

class FileReader{
public:
	FileReader(){}
	bool read(const string file_name, vector<string> &seq, vector<string> &seq_name);
};