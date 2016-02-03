CC = g++ $(CFLAGS)
CFLAGS = -g -O3

all: Holec

Holec: Holec.cpp input_function.cpp input_from_file.cpp Solve.cpp norms.cpp get_time.cpp functions.h
	g++ -o Holec Holec.cpp input_function.cpp input_from_file.cpp Solve.cpp norms.cpp get_time.cpp -pthread 

clean:
	rm -rf *.o Holec 
