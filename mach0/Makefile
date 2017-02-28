TARGET=mach0

.PHONY: utest vtest clean

$(TARGET): main.c
	gcc -o $@ $^ -lm

$(TARGET)_utest: main.c
	gcc -o $@ $^ -lm -DUTEST

$(TARGET)_vtest: main.c
	gcc -o $@ $^ -lm -DVTEST

utest: $(TARGET)_utest
	./$(TARGET)_utest

vtest: $(TARGET)_vtest
	./$(TARGET)_vtest

clean:
	-rm -f $(TARGET)
	-rm -f $(TARGET)_utest
	-rm -f $(TARGET)_vtest
