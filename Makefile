PLATFORM=nersc
include Makefile.in.$(PLATFORM)

.PHONY: clean all

all: centroid2d.x centroid2d++.x hello-mpi.x

%.x: %.c
	$(CC) $< -o $@

%.x: %.cc
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f centroid2d.x centroid2d++.x hello-mpi.x
