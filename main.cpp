#include "LinCapR.hpp"
#include "FileReader.hpp"

#include <iostream>
#include <fstream>
#include <cstring>

// Usage: ./LinCapR <input_file> <output_file> <beam_size> [options]
int main(int argc, char **argv){
	if(argc < 4){
		cout << "Usage: ./LinCapR <input_file> <output_file> <beam_size> [options]" << endl;
		return 1;
	}

	// args
	string input_file = argv[1];
	string output_file = argv[2];
	int beam_size = atoi(argv[3]);

	// get options
	bool output_energy = false;
	for(int i = 4; i < argc; i++){
		if(strcmp(argv[i], "-e") == 0){
			output_energy = true;
		}else{
			cout << "Error: invalid option: " << argv[i] << endl;
			return 1;
		}
	}

	// read fasta file
	FileReader fr;
	vector<string> seq, seq_name;
	fr.read(input_file, seq, seq_name);

	// check & clear output file
	ofstream ofs(output_file, ios::out | ios::trunc);
	if(!ofs){
		cout << "Error: cannot open output file: " << output_file << endl;
		return 1;
	}
	ofs.close();

	// run LinCapR
	const int s = seq.size();
	LinCapR lcr(beam_size);
	for(int i = 0; i < s; i++){
		// calc structural profile
		lcr.run(seq[i]);

		// output profile
		ofstream ofs(output_file, ios::out | ios::app);
		if(!ofs){
			cout << "Error: cannot open output file: " << output_file << endl;
			return 1;
		}
		lcr.output(ofs, seq_name[i]);
		ofs.close();

		if(output_energy) printf("G_ensemble: %.2lf\n", lcr.get_energy_ensemble());

		lcr.clear();
	}
}
