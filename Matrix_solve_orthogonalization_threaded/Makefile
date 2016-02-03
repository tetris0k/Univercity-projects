CC = g++ $(CFLAGS)
CFLAGS = -g -O3

all: doit

doit: main.cpp input_function.cpp input_from_file.cpp Solve.cpp norms.cpp get_time.cpp functions.h
	g++ -o doit main.cpp input_function.cpp input_from_file.cpp Solve.cpp norms.cpp get_time.cpp -pthread

clean:
	rm -rf *.o doit 
