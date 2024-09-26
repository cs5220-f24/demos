PLATFORM=nersc
include ../Makefile.in.$(PLATFORM)

.PHONY: clean all
EXE=mpi-demo.x

%.x: %.cc
	$(MPICXX) $(CXXFLAGS) $< -o $@

%.md: %.cc
	ldoc -o $@ $<

all: $(EXE)

clean:
	rm -f $(EXE)
