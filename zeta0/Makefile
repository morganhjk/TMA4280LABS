PROG=zeta0
COMP=gcc
FLAG=-lm -std=gnu11 -Wall -Wextra -Wpedantic

.PHONY: utest vtest clean

$(PROG): main.c
	$(COMP) -o $@ $^ $(FLAG)

$(PROG)_utest: main.c
	$(COMP) -o $@ $^ $(FLAG) -DUTEST

$(PROG)_vtest: main.c
	$(COMP) -o $@ $^ $(FLAG) -DVTEST

utest: $(PROG)_utest
	./$<

vtest: $(PROG)_vtest
	./$<

clean:
	-rm -f $(PROG)
	-rm -f $(PROG)_utest
	-rm -f $(PROG)_vtest
	-rm -f vtest.txt
