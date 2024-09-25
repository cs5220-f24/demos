PLATFORM=nersc
include ../Makefile.in.$(PLATFORM)

.PHONY: clean all
EXE=mpi-demo.x

%.x: %.cc
	$(MPICXX) $(CXXFLAGS) $< -o $@

all: $(EXE)

clean:
	rm -f $(EXE)
