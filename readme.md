### Development Notes 
- many changes happened in support_funcs.cpp, need to be carefull when check functions in UVZsolver.cpp since they call many kernels in support_funcs.cpp
#### Done:
- [x] Load all inputs from files
- [x] Tested Memcopy from host to device
-  [x] Rewrite makefile 

#### Todo:
- Get rid of global constant
 - [x] Make a struct that contain global constamt
 - [x] Copy the struct to device
 - [x] Using shared mem in each kernel to store the struct
 - [ ] Test kernels one by one after getting rid of constants
- Put things to gether
 - [ ] Call kernels one by one to check if things are in their right positions
