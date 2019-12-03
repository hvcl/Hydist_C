### Development Notes 
- many changes happened in support_funcs.cpp, need to be carefull when check functions in UVZsolver.cpp since they call many kernels in support_funcs.cpp
#### Done:
- [x] Get rid of global constant
 - [x] Load all inputs from files
 - [x] Tested Memcopy from host to device
 - [x] Rewrite makefile 
 - [x] Check kernels without constants in support_funcs.cpp
 - [x] Test kernels without constants in UVZsolver_multithreads.cpp

#### Todo:
- Put things to gether
 - [ ] Call kernels one by one to check if things are in their right positions
- Test final result
- Use new tridiag solver to test performance
