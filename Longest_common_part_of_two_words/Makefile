CC = g++ $(CFLAGS)
CFLAGS = -g -O3

all: start test opt

start: mp.cpp
	g++ -o start mp.cpp

test: tester.cpp
	g++ -o test tester.cpp

opt: optimized.cpp
	g++ -o opt optimized.cpp
clean:
	rm -rf *.o start test opt 
