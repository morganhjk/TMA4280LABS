TARGET=mach0

.PHONY: utest clean

$(TARGET): main.c
	gcc -o $@ $^ -lm

$(TARGET)_utest: main.c
	gcc -o $@ $^ -lm -DUTEST

utest: $(TARGET)_utest
	./$(TARGET)_utest

clean:
	-rm -f $(TARGET)
	-rm -f $(TARGET)_utest
