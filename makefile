CFLAGS = -I.
DEPS = engine.h 
OBJ = main.o engine.o support_funcs.o

all: $(OBJ)
	nvcc -arch=sm_60 $(OBJ) -o rep

%.o : %.cpp $(DEPS)
	nvcc -x cu -arch=sm_60 $(CFLAGS) -dc $< -o $@

clean:
	rm -f *.o 
