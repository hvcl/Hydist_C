CFLAGS = -I.
DEPS = engine.h support_funcs.h UVZSolver_multithread.h loader.h
OBJ = main.o engine.o support_funcs.o UVZSolver_multithread.o loader.cpp

all: $(OBJ)
	nvcc -arch=sm_61 $(OBJ) -o rep

%.o : %.cpp $(DEPS)
	nvcc -x cu -arch=sm_61 $(CFLAGS) -dc $< -o $@

clean:
	rm -f *.o 
