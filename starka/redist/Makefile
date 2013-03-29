#
# relies on 'IDS' definitions from parent makefile
#

IDS := BOOST_ID CODEMIN_ID SAMTOOLS_ID TABIX_ID

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(IDS), $(call assert-defined-indirect,$I))

id-to-dir = $($1)
DIRS := $(foreach I,$(IDS), $(call id-to-dir,$I))

.PHONY: $(IDS)

LINKS := boost samtools tabix

all: $(CODEMIN_ID) $(LINKS) 


$(IDS): $(DIRS)

safe-expand = \
mkdir -p $1.extract.tmp; \
$2 -dc < $1 | tar -xf - -C $1.extract.tmp; \
mv $1.extract.tmp/* .; \
rmdir $1.extract.tmp

maybe-expand = $(if $(wildcard $1),$(call safe-expand,$1,$2,$3),)

$(DIRS):
	$(call maybe-expand,$@.tar.bz2,bzip2,$@)
	$(call maybe-expand,$@.tar.gz,gzip,$@)


define install-boost
cd $(BOOST_ID); \
if ! test -f ./bjam; then \
  ./bootstrap.sh --with-libraries=program_options,test; \
fi; \
./bjam -a link=static; \
cd stage; \
rm -f include; \
ln -s .. include
endef


boost: $(BOOST_ID)
	$(install-boost) && cd $(CURDIR) && \
	ln -sf $< $@

samtools: $(SAMTOOLS_ID)
	$(MAKE) -C $(SAMTOOLS_ID) && \
	ln -sf $< $@

tabix: $(TABIX_ID)
	$(MAKE) -C $(TABIX_ID) && \
	ln -sf $< $@

clean:
	$(RM) -r $(DIRS) $(LINKS)

