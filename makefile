CC = nvcc -std=c++11
CFLAGS = -I. -x cu -arch=sm_61
DEPS = engine.h support_funcs.h UVZSolver_multithread.h loader.h sediment_transport.h
OBJ = main.o engine.o support_funcs.o UVZSolver_multithread.o loader.o sediment_transport.o

all: $(OBJ)
	$(CC) -arch=sm_61 $(OBJ) -o rep

%.o : %.cpp $(DEPS)
	$(CC) $(CFLAGS) -dc $< -o $@

clean:
	rm -f *.o 
