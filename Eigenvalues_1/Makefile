CC = g++ $(CFLAGS)
CFLAGS = -g -O3

all: start

start: main.cpp input_function.cpp input_from_file.cpp Solve.cpp norms.cpp functions.h
	g++ -g -o start main.cpp input_function.cpp input_from_file.cpp Solve.cpp norms.cpp

clean:
	rm -rf *.o start 
