PLATFORM=nersc
include ../Makefile.in.$(PLATFORM)

.PHONY: clean all
EXE=hello-mpi.x

%.x: %.cc
	$(MPICXX) $(CXXFLAGS) $< -o $@

all: $(EXE)

clean:
	rm -f $(EXE)
