
.PHONY: clean test

# TODO: find better system to apply all directories than foreach -- should be parallel, and catch an error in any one directory..

clean:
	$(foreach dir,$(SHARED_CXX_DIR), $(MAKE) -C $(dir) $@ &&) exit 0

test:
	$(foreach dir,$(SHARED_CXX_DIR), $(MAKE) -C $(dir) $@ &&) exit 0

