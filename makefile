CFLAGS = -I.
DEPS = engine.h support_funcs.h UVZSolver_multithread.h
OBJ = main.o engine.o support_funcs.o UVZSolver_multithread.o

all: $(OBJ)
	nvcc -arch=sm_60 $(OBJ) -o rep

%.o : %.cpp $(DEPS)
	nvcc -x cu -arch=sm_60 $(CFLAGS) -dc $< -o $@

clean:
	rm -f *.o 