CXX = g++
CXXFLAGS = -g -Wall -Wextra -O3 -lpthread

all: benchmark

benchmark: samples/benchmark.cpp src/pmc.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: clean

clean:
	rm -f benchmark
