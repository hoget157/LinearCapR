#include "FileReader.hpp"

#include <fstream>
#include <iostream>

// trim \r, \n in end of the line
inline void trim_end(string &s){
	while(!s.empty() && isspace(s.back())) s.pop_back();
}

// read sequences from file_name
bool FileReader::read(const string file_name, vector<string> &seq, vector<string> &seq_name){
	ifstream ifs(file_name);
	if(!ifs){
		cout << "Error: cannot open input file: " << file_name << endl;
		return false;
	}

	string line;
	while(getline(ifs, line)){
		trim_end(line);
		if(line.empty()) continue;
		if(line[0] == '>'){
			// new sequence
			seq_name.push_back(line.substr(1));
			seq.push_back("");
		}else{
			// append sequence
			seq.back() += line;
		}
	}

	ifs.close();
	return true;
}