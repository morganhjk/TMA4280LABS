TARGET=reduc

.PHONY: mpisum recursive clean

$(TARGET): main.c
	mpicc -o $@ $^ -lm -std=c99 -Wall -Wextra -Wpedantic

$(TARGET)_mpisum: main.c
	mpicc -o $@ $^ -lm -std=c99 -Wall -Wextra -DMPISUM

$(TARGET)_recursive: main.c
	mpicc -o $@ $^ -lm -std=c99 -Wall -Wextra -DRECURSIVE

mpisum: $(TARGET)_mpisum
	#./$(TARGET)_mpisum

recursive: $(TARGET)_recursive
	#./$(TARGET)_recursive

clean:
	-rm -f $(TARGET)
	-rm -f $(TARGET)_mpisum
	-rm -f $(TARGET)_recursive
