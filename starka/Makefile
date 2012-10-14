
SRC_DIR := $(CURDIR)/src
export BIN_DIR := $(CURDIR)/bin
export REDIST_DIR := $(CURDIR)/redist
export BOOST_ID := boost_1_44_0_subset
export CODEMIN_ID := CodeMin-1.0.2
export SAMTOOLS_ID := samtools-0.1.14_no_tview
export TABIX_ID := tabix-0.2.5_new_libz_order
export VCFTOOLS_ID := vcftools_0.1.6

all: install

build:
	$(MAKE) -C $(REDIST_DIR) && \
	$(MAKE) -C $(SRC_DIR)

install: build
	mkdir -p $(BIN_DIR) && $(MAKE) -C $(SRC_DIR) $@ 

clean: srcclean
	$(MAKE) -C $(REDIST_DIR) $@ 

# developer targets:
#

# clean only src directories but leave redist built:
srcclean:
	$(MAKE) -C $(SRC_DIR) clean
	$(RM) -rf $(BIN_DIR)

etags:
	$(MAKE) -C $(SRC_DIR) $@

