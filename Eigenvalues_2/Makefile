CC = g++ $(CFLAGS)
CFLAGS = -g -O3

all: doit

doit: main.cpp input_function.cpp input_from_file.cpp Solve.cpp functions.h
	g++ -o doit main.cpp input_function.cpp input_from_file.cpp Solve.cpp

clean:
	rm -rf *.o doit 
