CC = g++ $(CFLAGS)
CFLAGS = -g -O3

all: doit test new

doit: matprak.cpp
	g++ -o doit matprak.cpp

test: test.cpp
	g++ -o test test.cpp

new: noviy.cpp
	g++ -o new noviy.cpp

clean:
	rm -rf *.o doit test new
