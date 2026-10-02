CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -pthread -Wall -Wextra

primes: primes.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	del /Q primes

.PHONY: clean
