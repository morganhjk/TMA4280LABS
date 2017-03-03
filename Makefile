MSC=$(MAKE) -s -C

.PHONY: zeta0 zeta1 mach0 mach1 utest vtest clean

zeta0:
	@$(MSC) zeta0

zeta1:
	@$(MSC) zeta1

mach0:
	@$(MSC) mach0

mach1:
	@$(MSC) mach1

utest:
	@$(MSC) zeta0 utest
	@$(MSC) mach0 utest

vtest:
	@$(MSC) zeta0 vtest
	@$(MSC) mach0 vtest
	@$(MSC) zeta1 vtest
	@$(MSC) mach1 vtest

clean:
	@$(MSC) zeta0 clean
	@$(MSC) zeta1 clean
	@$(MSC) mach0 clean
	@$(MSC) mach1 clean
