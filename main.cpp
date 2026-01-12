#include "LinCapR.hpp"
#include "FileReader.hpp"
#include "energy_raccess.hpp"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

// Usage: ./LinCapR <input_file> <output_file> <beam_size> [options]
int main(int argc, char **argv){
	if(argc < 4){
		cout << "Usage: ./LinCapR <input_file> <output_file> <beam_size> [options]" << endl;
		cout << "Options:" << endl;
		cout << "  -e                 Output ensemble energy" << endl;
		cout << "  --logz             Output logZ (alpha_O[last])" << endl;
		cout << "  --compare-unpaired Compare 1 - p(stem) with Raccess p([i:i], unpaired)" << endl;
		cout << "  --max-span <int>   Raccess max span for compare-unpaired (default: seq length)" << endl;
		cout << "  --energy <model>   Energy model: turner2004 (default) or turner1999" << endl;
		cout << "  --engine <name>    Energy engine: lincapr (default) or raccess" << endl;
		return 1;
	}

	// args
	string input_file = argv[1];
	string output_file = argv[2];
	int beam_size = atoi(argv[3]);

	// get options
	bool output_energy = false;
	bool output_logz = false;
	bool compare_unpaired = false;
	int compare_max_span = 0;
	energy::Model energy_model = energy::Model::Turner2004;
	LinCapR::EnergyEngine engine = LinCapR::EnergyEngine::LinearCapR;
	for(int i = 4; i < argc; i++){
		if(strcmp(argv[i], "-e") == 0){
			output_energy = true;
		}else if(strcmp(argv[i], "--logz") == 0){
			output_logz = true;
		}else if(strcmp(argv[i], "--compare-unpaired") == 0){
			compare_unpaired = true;
		}else if(strcmp(argv[i], "--max-span") == 0){
			if(i + 1 >= argc){
				cout << "Error: --max-span requires an integer argument" << endl;
				return 1;
			}
			compare_max_span = atoi(argv[++i]);
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
		}else if(strcmp(argv[i], "--engine") == 0){
			if(i + 1 >= argc){
				cout << "Error: --engine requires an argument (lincapr or raccess)" << endl;
				return 1;
			}
			const char *choice = argv[++i];
			if(strcmp(choice, "lincapr") == 0){
				engine = LinCapR::EnergyEngine::LinearCapR;
			}else if(strcmp(choice, "raccess") == 0){
				engine = LinCapR::EnergyEngine::Raccess;
			}else{
				cout << "Error: invalid engine: " << choice << endl;
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
		}else if(strncmp(argv[i], "--engine=", 9) == 0){
			const char *choice = argv[i] + 9;
			if(strcmp(choice, "lincapr") == 0){
				engine = LinCapR::EnergyEngine::LinearCapR;
			}else if(strcmp(choice, "raccess") == 0){
				engine = LinCapR::EnergyEngine::Raccess;
			}else{
				cout << "Error: invalid engine: " << choice << endl;
				return 1;
			}
		}else if(strncmp(argv[i], "--max-span=", 11) == 0){
			compare_max_span = atoi(argv[i] + 11);
		}else{
			cout << "Error: invalid option: " << argv[i] << endl;
			return 1;
		}
	}
	if(compare_unpaired && engine != LinCapR::EnergyEngine::Raccess){
		cout << "Error: --compare-unpaired requires --engine raccess" << endl;
		return 1;
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
	LinCapR lcr(beam_size, energy_model, engine);
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
