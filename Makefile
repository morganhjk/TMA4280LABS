MSC=$(MAKE) -s -C

.PHONY: zeta0 zeta1 zeta2 mach0 mach1 mach2 reduc utest vtest clean

zeta0:
	@$(MSC) $@

zeta1:
	@$(MSC) $@

zeta2:
	@$(MSC) $@

mach0:
	@$(MSC) $@

mach1:
	@$(MSC) $@

mach2:
	@$(MSC) $@

reduc:
	@$(MSC) $@

utest:
	@$(MSC) zeta0 $@
	@$(MSC) mach0 $@

vtest:
	@$(MSC) zeta0 $@
	@$(MSC) zeta1 $@
	@$(MSC) zeta2 $@
	@$(MSC) mach0 $@
	@$(MSC) mach1 $@
	@$(MSC) mach2 $@

clean:
	@$(MSC) zeta0 $@
	@$(MSC) zeta1 $@
	@$(MSC) zeta2 $@
	@$(MSC) mach0 $@
	@$(MSC) mach1 $@
	@$(MSC) mach2 $@
	@$(MSC) reduc $@
