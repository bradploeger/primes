CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -pthread -Wall -Wextra

primes: primes.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f primes

.PHONY: clean
