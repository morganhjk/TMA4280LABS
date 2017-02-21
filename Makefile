.PHONY: zeta0 zeta1 mach0 mach1 utest clean

zeta0:
	@$(MAKE) -s -C zeta0

zeta1:
	@$(MAKE) -s -C zeta1

mach0:
	@$(MAKE) -s -C mach0

mach1:
	@$(MAKE) -s -C mach1

utest:
	@$(MAKE) -s -C zeta0 utest
	@$(MAKE) -s -C mach0 utest

clean:
	@$(MAKE) -s -C zeta0 clean
	@$(MAKE) -s -C zeta1 clean
	@$(MAKE) -s -C mach0 clean
	@$(MAKE) -s -C mach1 clean
