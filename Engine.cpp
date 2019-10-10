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


pair<int, int> load_inputs(string dir,
				vector<DOUBLE> &h,
				vector<DOUBLE> &hsnham,
				vector<DOUBLE> &VISCOINDX, 
				vector<DOUBLE> &bc_up,
				vector<DOUBLE> &bc_right,
				vector<DOUBLE> &bc_left, 
				vector<DOUBLE> &bc_down, 
				vector<DOUBLE> &CC_u, 
				vector<DOUBLE> &CC_d,
				vector<DOUBLE> &CC_l,
				vector<DOUBLE> &CC_r,
				vector<int> &bienQ)
{
	pair<int, int > size;
	size = load_file<DOUBLE> (dir + "bandodosau.txt", h);
	// cout << size.first << " " << size.second << endl;
	load_file<DOUBLE> (dir + "hsnham.txt", hsnham);

	cout << size.first << " " << size.second << endl;
	load_file<DOUBLE> (dir + "hsnhotroiA.txt", VISCOINDX);

	// load_file<DOUBLE> (dir + "bientren.txt", bc_up);
	load_file<DOUBLE> (dir + "bienduoi.txt", bc_down);
	// load_file<DOUBLE> (dir + "bientrai.txt", bc_left);
	load_file<DOUBLE> (dir + "bienphai.txt", bc_right);
	// load_file<DOUBLE> (dir + "FSbientren.txt", CC_u);
	// load_file<DOUBLE> (dir + "FSbienduoi.txt", CC_d);
	// load_file<DOUBLE> (dir + "FSbientrai.txt", CC_l);
	// load_file<DOUBLE> (dir + "FSbienphai.txt", CC_r);
	load_file<int> (dir + "boundary_type.txt", bienQ);
	return size;
}

void load_initial_condition(string dir, 
							vector<DOUBLE> &u, 
							vector<DOUBLE> &v,
							vector<DOUBLE> &z,
							vector<DOUBLE> &FS,
							vector<int> &khouot)
{

	load_file<DOUBLE> (dir + "u.txt", u);
	load_file<DOUBLE> (dir + "v.txt", v);
	load_file<DOUBLE> (dir + "z.txt", z);
	load_file<DOUBLE> (dir + "FS.txt", FS);
	load_file<int> (dir + "khouot.txt", khouot);
}