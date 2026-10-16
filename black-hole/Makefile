CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra

blackhole: blackhole.cpp
	$(CXX) $(CXXFLAGS) -pthread -o $@ $<

test: blackhole
	./blackhole --test

clean:
	rm -f blackhole

.PHONY: test clean
