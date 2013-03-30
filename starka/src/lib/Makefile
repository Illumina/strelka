
.PHONY: clean test

clean:
	for dir in $(SHARED_CXX_DIR); do \
		$(MAKE) -C $$dir $@; \
	done

test:
	for dir in $(SHARED_CXX_DIR); do \
		$(MAKE) -C $$dir $@; \
	done

