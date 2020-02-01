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

#### Running Instruction:  
 1. Fill in correct coefficients for input mesh in Load_coeffs() function in loader.cpp files
 2. Copy input files to Inputs/ directory with corresponding names: 
  - bandodosau.txt : depth map of to-be compute mesh data
  - hsnhotroiA.txt: visco coefficient map
  - hsnham.txt: roughness index map
  - Fw_map.txt: Fw coefficient map
  - FS_bientren.txt: FS boundary condition for upper boundary
  - FS_bienduoi.txt: FS boundary condition for lower boundary
  - FS_bientrai.txt: FS boundary condition for left side boundary
  - FS_bienphai.txt: FS boundary condition for right side boundary
  - bientren.txt: boundary condition for upper boundary
  - bienduoi.txt: boundary condition for lower boundary
  - bientrai.txt: boundary condition for left side boundary
  - bienphai.txt: boundary condition for right side boundary
  - boundary_type.txt: boundary type of each side (solid type or water type) 
 3. Run standardize_inp.py to transpose inputs into Python-style coordinates
 4. run command:  ```console ~$ make ``` to compile changes
 5. input arguments in args.txt with following options: 
  - -f Directory that contains input files 
    -h Number of running hour
    -m number of running mins
    -s number of running seconds
    -i time interval in which output is saved
    -b time, in second, at which Sediment transport is started to be computed
 6. run command: ``` console ./rep `cat args.txt` ``` to run rep file with the specified arguments in args.txt