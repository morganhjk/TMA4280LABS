TARGET=mach1

.PHONY: clean

$(TARGET): main.c
	mpicc -o $@ $^ -lm -std=c99 -Wall -Wextra -Wpedantic

clean:
	-rm -f $(TARGET)
