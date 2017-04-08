PROG0=mach1
PROG1=$(PROG0)_vtest

CC=mpicc
CFLAGS=-lm -std=gnu11 -Wall -Wextra -Wpedantic

.PHONY: vtest clean

$(PROG0): main.c
	$(CC) -o $@ $^ $(CFLAGS)

$(PROG1): main.c
	$(CC) -o $@ $^ $(CFLAGS) -DVTEST

vtest: $(PROG1)
	mpirun -np 2 $<
	mpirun -np 4 $<
	mpirun -np 8 $<

clean:
	-rm -f $(PROG0)
	-rm -f $(PROG1)
	-rm -f vtest.txt
