CC = mpicxx $(CFLAGS)
CFLAGS = -g -O3

all: doit

doit: main.cpp input_function.cpp input_from_file.cpp Solve.cpp norms.cpp functions.h
	mpicxx -o doit main.cpp input_function.cpp input_from_file.cpp Solve.cpp norms.cpp
clean:
	rm -rf *.o doit 
