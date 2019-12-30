### Development Notes 
- many changes happened in support_funcs.cpp, need to be carefull when check functions in UVZsolver.cpp since they call many kernels in support_funcs.cpp
#### Done:
- Get rid of global constant
 - [x] Load all inputs from files
 - [x] Tested Memcopy from host to device
 - [x] Rewrite makefile 
 - [x] Check kernels without constants in support_funcs.cpp
 - [x] Test kernels without constants in UVZsolver_multithreads.cpp
 - [x] change sediment transport kernels to work with non_constant versions
- Put things to gether
 - [x] Call kernels one by one to check if things are in their right positions
 - [x] Test final result for hydraulic module
 - [x] test sediment transport functions
 - [x] test the program with sediment transport

#### Todo:

- [ ] Use new tridiag solver to test performance

