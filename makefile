CC = nvcc -std=c++11
FLAGS = -I. -x cu -arch=sm_61 --ptxas-options=-v -Xptxas -dlcm=ca -O3


DEPS = engine.h support_funcs.h UVZSolver_multithread.h loader.h sediment_transport.h gtsv.h
OBJ = main.o engine.o support_funcs.o UVZSolver_multithread.o loader.o sediment_transport.o gtsv.o

all: $(OBJ)
	$(CC) -arch=sm_61 -lcusparse $(OBJ) -o rep

%.o : %.cpp $(DEPS)
	$(CC) $(FLAGS) -dc $< -o $@

clean:
	rm -f *.o 
