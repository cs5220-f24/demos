PLATFORM=nersc
include ../Makefile.in.$(PLATFORM)

.PHONY: clean all
EXE=ping.x

%.x: %.c
	$(MPICC) $(CFLAGS) $< -o $@

all: $(EXE)

clean:
	rm -f $(EXE)
