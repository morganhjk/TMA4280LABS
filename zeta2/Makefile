PROG0=zeta2
PROG1=$(PROG0)_vtest
PROG2=$(PROG0)_lille

CC=mpicc
CFLAGS=-lm -fopenmp -std=gnu11 -Wall -Wextra -Wpedantic

.PHONY: vtest clean

$(PROG0): main.c
	$(CC) -o $@ $^ $(CFLAGS)

$(PROG1): main.c
	$(CC) -o $@ $^ $(CFLAGS) -DVTEST

$(PROG2): lille.c
	$(CC) -o $@ $^ $(CFLAGS) -D_GNU_SOURCE

vtest: $(PROG1)
	mpirun -np 2 $<
	mpirun -np 4 $<
	mpirun -np 8 $<

clean:
	-rm -f $(PROG0)
	-rm -f $(PROG1)
	-rm -f $(PROG2)
	-rm -f vtest.txt
	-rm -f lille-*.txt
