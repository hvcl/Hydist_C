#pragma once
#include <map>
#include <string>
#include <algorithm>
#include "cuda_runtime.h"
using namespace std;

#define DOUBLE double

struct pointer_struct {
	double* ptr;
	int size;
	pointer_struct() {
		ptr = NULL;
		size = 0;
	}
	pointer_struct(double* p, int s) {
		ptr = p;
		size = s;
	}
	
};
class Pointer
{
public:
	Pointer(map<string, pointer_struct> &ptr_list, int M, int N);
	~Pointer();
private:
	map<string, pointer_struct> host_pointers;
	map<string, pointer_struct> device_pointers;
	void to_device();
	void extract(pointer_struct* arg_list[], int size);
	double* get_device_pointer(string key);
	string find_value(double* ptr);
};

