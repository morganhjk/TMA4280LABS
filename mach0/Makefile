TARGET=mach0

.PHONY: utest clean

$(TARGET): main.c
	gcc -o $@ $^ -lm

utest: $(TARGET) test-unit.sh
	./test-unit.sh

clean:
	-rm -f $(TARGET)
