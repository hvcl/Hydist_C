CFLAGS = -I.
DEPS = Engine.h 
OBJ = main.o Engine.o support_funcs.o

all: $(OBJ)
	nvcc -arch=sm_60 $(OBJ) -o rep

%.o : %.cpp $(DEPS)
	nvcc -x cu -arch=sm_60 -I. -dc $< -o $@

clean:
	rm -f *.o 
