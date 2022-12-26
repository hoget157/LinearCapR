#include "LinearCapR.hpp"
#include "FileReader.hpp"
#include "utils.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

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

	// run LinearCapR
	const int s = seq.size();
	LinearCapR lcr(beam_size);
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

		lcr.clear();
	}
}

// int main(){
// 	LinearCapR lcr(100);

// 	int n;
// 	cin >> n;
// 	vector<Float> vec(n);
// 	for(int i = 0;i < n;i++) cin >> vec[i];
// 	for(int i = 0;i < n;i++){
// 		vector<Float> vec2 = vec;
// 		cout << lcr.quickselect(vec2, 0, n, i + 1) << endl;
// 	}
// }