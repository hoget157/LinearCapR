#include <iostream>
#include <string>
#include <vector>

#include "LinearCapR.hpp"
#include "FileReader.hpp"
#include "utils.hpp"

using namespace std;

// Usage: ./LinearCapR <input_file> <output_file> <beam_size>
int main(int argc, char **argv){
	if(argc != 4){
		cout << "Usage: ./LinearCapR <input_file> <output_file> <beam_size>" << endl;
		return 1;
	}

	// args
	string input_file = argv[1];
	string output_file = argv[2];
	int beam_size = atoi(argv[3]);

	// read fasta file
	FileReader fr;
	vector<string> seq, seq_name;
	fr.read(input_file, seq, seq_name);

	cout << seq_name << endl;
	cout << seq << endl;

	// run LinearCapR
// 	const int s = seq.size();
// 	LinearCapR lcr(beam_size);
// 	for(int i = 0; i < s; i++){
// 		// calc structural profile
// 		lcr.run(seq[i]);

// 		// output profile
// 		cout << ">" << seq_name[i] << endl;
		
// ofstream ofs(out_file);
// 	if(!ofs){
// 		cout << "cannot open output file" << endl;
// 		exit(1);
// 	}
// 	ofs << ">genome" << endl;
// 	ofs << genome << endl;
// 	}
}