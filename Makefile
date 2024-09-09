CC=cc
CXX=c++
CXXFLAGS=-std=gnu++20

.PHONY: clean all

all: centroid2d.x centroid2d++.x

centroid2d.x: centroid2d.c
	$(CC) $< -o $@

centroid2d++.x: centroid2d.cc
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f centroid2d.x centroid2d++.x
