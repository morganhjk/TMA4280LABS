PROG0=reduc
PROG1=$(PROG0)_mpisum
PROG2=$(PROG0)_recursive

CC=mpicc
CFLAGS=-lm -std=gnu11 -Wall -Wextra -Wpedantic

.PHONY: all test clean

all: $(PROG0) $(PROG1) $(PROG2)

$(PROG0): main.c
	$(CC) -o $@ $^ $(CFLAGS)

$(PROG1): main.c
	$(CC) -o $@ $^ $(CFLAGS) -DMPISUM

$(PROG2): main.c
	$(CC) -o $@ $^ $(CFLAGS) -DRECURSIVE

test: $(PROG0) $(PROG1) $(PROG2)
	mpirun -np 2 $(PROG0) $(N)
	mpirun -np 4 $(PROG0) $(N)
	mpirun -np 8 $(PROG0) $(N)
	mpirun -np 2 $(PROG1) $(N)
	mpirun -np 4 $(PROG1) $(N)
	mpirun -np 8 $(PROG1) $(N)
	mpirun -np 2 $(PROG2) $(N)
	mpirun -np 4 $(PROG2) $(N)
	mpirun -np 8 $(PROG2) $(N)

clean:
	-rm -f $(PROG0)
	-rm -f $(PROG1)
	-rm -f $(PROG2)
