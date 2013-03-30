
.PHONY: clean test

clean:
	$(foreach dir,$(SHARED_CXX_DIR), $(MAKE) -C $(dir) $@ &&) exit 0

test:
	$(foreach dir,$(SHARED_CXX_DIR), $(MAKE) -C $(dir) $@ &&) exit 0

