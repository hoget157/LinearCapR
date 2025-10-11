#include "LinCapR.hpp"
#include "FileReader.hpp"

#include <iostream>
#include <fstream>
#include <cstring>

// Usage: ./LinCapR <input_file> <output_file> <beam_size> [options]
int main(int argc, char **argv){
	if(argc < 4){
		cout << "Usage: ./LinCapR <input_file> <output_file> <beam_size> [options]" << endl;
		cout << "Options:" << endl;
		cout << "  -e                 Output ensemble energy" << endl;
		cout << "  --energy <model>   Energy model: turner2004 (default) or turner1999" << endl;
		return 1;
	}

	// args
	string input_file = argv[1];
	string output_file = argv[2];
	int beam_size = atoi(argv[3]);

	// get options
	bool output_energy = false;
	energy::Model energy_model = energy::Model::Turner2004;
	for(int i = 4; i < argc; i++){
		if(strcmp(argv[i], "-e") == 0){
			output_energy = true;
		}else if(strcmp(argv[i], "--energy") == 0){
			if(i + 1 >= argc){
				cout << "Error: --energy requires an argument (turner2004 or turner1999)" << endl;
				return 1;
			}
			const char *choice = argv[++i];
			if(strcmp(choice, "turner2004") == 0){
				energy_model = energy::Model::Turner2004;
			}else if(strcmp(choice, "turner1999") == 0){
				energy_model = energy::Model::Turner1999;
			}else{
				cout << "Error: invalid energy model: " << choice << endl;
				return 1;
			}
		}else if(strncmp(argv[i], "--energy=", 9) == 0){
			const char *choice = argv[i] + 9;
			if(strcmp(choice, "turner2004") == 0){
				energy_model = energy::Model::Turner2004;
			}else if(strcmp(choice, "turner1999") == 0){
				energy_model = energy::Model::Turner1999;
			}else{
				cout << "Error: invalid energy model: " << choice << endl;
				return 1;
			}
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
	LinCapR lcr(beam_size, energy_model);
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
