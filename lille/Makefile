PROG=lille

CC=mpicc
CFLAGS=-lm -fopenmp -std=gnu11 -Wall -Wextra -Wpedantic -D_GNU_SOURCE

.PHONY: clean

$(PROG): main.c
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	-rm -f $(PROG)
	-rm -f $(PROG)-*.txt
