PROG=zeta1
COMP=mpicc
FLAG=-lm -std=gnu11 -Wall -Wextra -Wpedantic

.PHONY: vtest clean

$(PROG): main.c
	$(COMP) -o $@ $^ $(FLAG)

$(PROG)_vtest: main.c
	$(COMP) -o $@ $^ $(FLAG) -DVTEST

vtest: $(PROG)_vtest
	mpirun -np 2 $<
	mpirun -np 4 $<
	mpirun -np 8 $<

clean:
	-rm -f $(PROG)
	-rm -f $(PROG)_vtest
	-rm -f vtest.txt
