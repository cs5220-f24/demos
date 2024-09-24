PLATFORM=nersc
include ../Makefile.in.$(PLATFORM)

EXE=centroid2d.x centroid2d++.x
.PHONY: clean all

%.x: %.c
	$(CC) $(CFLAGS) $< -o $@

%.x: %.cc
	$(CXX) $(CXXFLAGS) $< -o $@ 

all: $(EXE)

clean:
	rm -f $(EXE)
