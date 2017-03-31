PROG0=mach0
PROG1=$(PROG0)_utest
PROG2=$(PROG0)_vtest

CC=gcc
CFLAGS=-lm -std=gnu11 -Wall -Wextra -Wpedantic

.PHONY: utest vtest clean

$(PROG0): main.c
	$(CC) -o $@ $^ $(CFLAGS)

$(PROG1): main.c
	$(CC) -o $@ $^ $(CFLAGS) -DUTEST

$(PROG2): main.c
	$(CC) -o $@ $^ $(CFLAGS) -DVTEST

utest: $(PROG1)
	./$<

vtest: $(PROG2)
	./$<

clean:
	-rm -f $(PROG0)
	-rm -f $(PROG1)
	-rm -f $(PROG2)
	-rm -f vtest.txt
