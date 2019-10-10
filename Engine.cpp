#include "Engine.h"

template <typename T>
pair<int, int> load_file(string filename, vector<T> arr){

	ifstream ifs;
	ifs.open(filename);
	if (!ifs) {
		cout << "cannot open" <<  filename << endl;
	}
	string line;
	int M, N;
	N = 0;
	while (getline(ifs, line)){
		N++;
		stringstream linestream(line);
		T x;
		M = 0;
		while (linestream >> x){
			M++;
			arr.push_back(x);
		}
		arr.push_back(arr[arr.size() - 1]); // pad edge
	}
	return pair<int, int>(N, M);
}

template pair<int, int> load_file<double>(string filename, std::vector<double> arr);
template pair<int, int> load_file<float>(string filename, std::vector<float> arr);


