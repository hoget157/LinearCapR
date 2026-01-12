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
		cout << "  --logz-raccess     Output Raccess logZ (ProbModel partition coeff)" << endl;
		cout << "  --compare-unpaired Compare 1 - p(stem) with Raccess p([i:i], unpaired)" << endl;
		cout << "  --max-span <int>   Raccess max span for compare-unpaired (default: seq length)" << endl;
		cout << "  --energy <model>   Energy model: turner2004 (default) or turner1999" << endl;
		cout << "  --engine <name>    Energy engine: lincapr (default) or raccess" << endl;
		cout << "  --debug-local i,j  Compare closed-wrapper vs direct Raccess energies at (i,j) (0-origin)" << endl;
		cout << "  --debug-loop p,q   Optional inner pair for loop energy check (0-origin)" << endl;
		cout << "  --debug-stem i     Print top paired positions contributing to Stem[i]" << endl;
		cout << "  --debug-top n      Number of top pairs to print (default: 10)" << endl;
		cout << "  --debug-pair i,j   Print alpha/beta/prob for pair (i,j) (0-origin)" << endl;
		return 1;
	}

	// args
	string input_file = argv[1];
	string output_file = argv[2];
	int beam_size = atoi(argv[3]);

	// get options
	bool output_energy = false;
	bool output_logz = false;
	bool output_logz_raccess = false;
	bool compare_unpaired = false;
	int compare_max_span = 0;
	bool debug_local = false;
	bool debug_loop = false;
	int debug_i = -1;
	int debug_j = -1;
	int debug_p = -1;
	int debug_q = -1;
	bool debug_stem = false;
	int debug_stem_i = -1;
	int debug_top = 10;
	bool debug_pair = false;
	int debug_pair_i = -1;
	int debug_pair_j = -1;
	energy::Model energy_model = energy::Model::Turner2004;
	LinCapR::EnergyEngine engine = LinCapR::EnergyEngine::LinearCapR;
	for(int i = 4; i < argc; i++){
		if(strcmp(argv[i], "-e") == 0){
			output_energy = true;
		}else if(strcmp(argv[i], "--logz") == 0){
			output_logz = true;
		}else if(strcmp(argv[i], "--logz-raccess") == 0){
			output_logz_raccess = true;
		}else if(strcmp(argv[i], "--compare-unpaired") == 0){
			compare_unpaired = true;
		}else if(strcmp(argv[i], "--max-span") == 0){
			if(i + 1 >= argc){
				cout << "Error: --max-span requires an integer argument" << endl;
				return 1;
			}
			compare_max_span = atoi(argv[++i]);
		}else if(strcmp(argv[i], "--debug-local") == 0){
			if(i + 1 >= argc){
				cout << "Error: --debug-local requires i,j" << endl;
				return 1;
			}
			const string arg = argv[++i];
			const size_t comma = arg.find(',');
			if(comma == string::npos){
				cout << "Error: --debug-local expects i,j (0-origin)" << endl;
				return 1;
			}
			debug_i = stoi(arg.substr(0, comma));
			debug_j = stoi(arg.substr(comma + 1));
			debug_local = true;
		}else if(strcmp(argv[i], "--debug-loop") == 0){
			if(i + 1 >= argc){
				cout << "Error: --debug-loop requires p,q" << endl;
				return 1;
			}
			const string arg = argv[++i];
			const size_t comma = arg.find(',');
			if(comma == string::npos){
				cout << "Error: --debug-loop expects p,q (0-origin)" << endl;
				return 1;
			}
			debug_p = stoi(arg.substr(0, comma));
			debug_q = stoi(arg.substr(comma + 1));
			debug_loop = true;
		}else if(strcmp(argv[i], "--debug-stem") == 0){
			if(i + 1 >= argc){
				cout << "Error: --debug-stem requires i" << endl;
				return 1;
			}
			debug_stem_i = atoi(argv[++i]);
			debug_stem = true;
		}else if(strcmp(argv[i], "--debug-top") == 0){
			if(i + 1 >= argc){
				cout << "Error: --debug-top requires n" << endl;
				return 1;
			}
			debug_top = atoi(argv[++i]);
		}else if(strcmp(argv[i], "--debug-pair") == 0){
			if(i + 1 >= argc){
				cout << "Error: --debug-pair requires i,j" << endl;
				return 1;
			}
			const string arg = argv[++i];
			const size_t comma = arg.find(',');
			if(comma == string::npos){
				cout << "Error: --debug-pair expects i,j (0-origin)" << endl;
				return 1;
			}
			debug_pair_i = stoi(arg.substr(0, comma));
			debug_pair_j = stoi(arg.substr(comma + 1));
			debug_pair = true;
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
		}else if(strncmp(argv[i], "--debug-local=", 14) == 0){
			const string arg = argv[i] + 14;
			const size_t comma = arg.find(',');
			if(comma == string::npos){
				cout << "Error: --debug-local expects i,j (0-origin)" << endl;
				return 1;
			}
			debug_i = stoi(arg.substr(0, comma));
			debug_j = stoi(arg.substr(comma + 1));
			debug_local = true;
		}else if(strncmp(argv[i], "--debug-loop=", 13) == 0){
			const string arg = argv[i] + 13;
			const size_t comma = arg.find(',');
			if(comma == string::npos){
				cout << "Error: --debug-loop expects p,q (0-origin)" << endl;
				return 1;
			}
			debug_p = stoi(arg.substr(0, comma));
			debug_q = stoi(arg.substr(comma + 1));
			debug_loop = true;
		}else if(strncmp(argv[i], "--debug-stem=", 13) == 0){
			debug_stem_i = atoi(argv[i] + 13);
			debug_stem = true;
		}else if(strncmp(argv[i], "--debug-top=", 12) == 0){
			debug_top = atoi(argv[i] + 12);
		}else if(strncmp(argv[i], "--debug-pair=", 13) == 0){
			const string arg = argv[i] + 13;
			const size_t comma = arg.find(',');
			if(comma == string::npos){
				cout << "Error: --debug-pair expects i,j (0-origin)" << endl;
				return 1;
			}
			debug_pair_i = stoi(arg.substr(0, comma));
			debug_pair_j = stoi(arg.substr(comma + 1));
			debug_pair = true;
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
		if(output_logz) printf("logZ: %.10f\n", lcr.get_logZ());
		if(output_logz_raccess){
			const int max_span = (compare_max_span > 0 ? compare_max_span : static_cast<int>(seq[i].size() + 1));
			const double logz_raccess = lcr::compute_raccess_logz(seq[i], max_span);
			printf("raccess_logZ: %.10f\n", logz_raccess);
			if(output_logz){
				printf("logZ_diff: %.10g\n", lcr.get_logZ() - logz_raccess);
			}
		}
		if(compare_unpaired){
			const int max_span = (compare_max_span > 0 ? compare_max_span : static_cast<int>(seq[i].size() + 1));
			const auto unpaired = lcr::compute_raccess_unpaired_1(seq[i], max_span);
			const auto& p_stem = lcr.get_prob_stem();
			double max_diff = -1.0;
			int max_i = -1;
			double max_unp = 0.0;
			double max_from_stem = 0.0;
			for(size_t k = 0; k < p_stem.size(); k++){
				const double p_unp_from_stem = 1.0 - p_stem[k];
				const double diff = fabs(p_unp_from_stem - unpaired[k]);
				if(diff > max_diff){
					max_diff = diff;
					max_i = static_cast<int>(k);
					max_unp = unpaired[k];
					max_from_stem = p_unp_from_stem;
				}
			}
			printf("unpaired_check: max_abs_diff=%.6g at i=%d (raccess_unpaired=%.6g, 1-p_stem=%.6g)\n",
			       max_diff, max_i, max_unp, max_from_stem);
		}

		if(debug_local){
			if(engine != LinCapR::EnergyEngine::Raccess){
				cout << "Warning: --debug-local uses the Raccess energy model regardless of --engine" << endl;
			}
			lcr::debug_raccess_local(seq[i], debug_i, debug_j, debug_loop, debug_p, debug_q);
		}
		if(debug_stem){
			lcr.debug_stem_pairs(debug_stem_i, debug_top);
		}
		if(debug_pair){
			lcr.debug_pair(debug_pair_i, debug_pair_j);
		}

		lcr.clear();
	}
}
