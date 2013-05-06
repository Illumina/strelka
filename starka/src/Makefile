
# Assumes values have been set by parent makefile:
#
# BOOST_ROOT : currently tested on Boost v1.44
#
# REDIST_DIR : location of other dependencies (besides boost)
#
# SAMTOOLS_ID : Uses the samtools library for bam parsing. This is included in redist/ and should not need to be set.
# TABIX_ID :
# CODEMIN_ID : 
#

IDS := REDIST_DIR BOOST_ID CODEMIN_ID SAMTOOLS_ID TABIX_ID

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(IDS), $(call assert-defined-indirect,$I))

file-exists = $(wildcard $1)
check-exists-indirect = $(if $(call file-exists,$($1)),,$(error Variable '$1' refers to a non-existent path '$($1)'))
$(foreach I,REDIST_DIR, $(call check-exists-indirect,$I))

SAMTOOLS_DIR := $(REDIST_DIR)/$(SAMTOOLS_ID)
TABIX_DIR := $(REDIST_DIR)/$(TABIX_ID)
CODEMIN_INCLUDE_DIR := $(REDIST_DIR)/$(CODEMIN_ID)/include

CXX := g++
BASE_WARNFLAGS := -Wall -Wextra -Wshadow -Wunused -g
EXTRA_WARNFLAGS := -Wpointer-arith -Winit-self -Wredundant-decls # -Wsign-conversion -Weverything 
WARNFLAGS = $(BASE_WARNFLAGS) $(EXTRA_WARNFLAGS)

STRICT_WARNFLAGS := -Werror
ifeq ($(STRICT_COMPILE),1)
    WARNFLAGS += $(STRICT_WARNFLAGS)
endif

OPTFLAGS := -O2
DEBUGFLAGS := -ggdb -O0
XFLAGS := $(WARNFLAGS) $(OPTFLAGS)

BOOST_REDIST_ROOT := $(REDIST_DIR)/$(BOOST_ID)/stage
BOOST_INCLUDE_DIR := $(BOOST_REDIST_ROOT)/include
BOOST_LIB_DIR := $(BOOST_REDIST_ROOT)/lib

STARKA_INCLUDE_DIR:=$(CURDIR)/lib

# use isystem so that we can run starka with stringent warnings and skip
# all the failed cases in boost
export CXXFLAGS := $(XFLAGS) -isystem $(BOOST_INCLUDE_DIR) -isystem $(SAMTOOLS_DIR) -isystem $(TABIX_DIR) -isystem $(CODEMIN_INCLUDE_DIR) -I$(STARKA_INCLUDE_DIR)
export LDFLAGS := -L$(BOOST_LIB_DIR) -L$(SAMTOOLS_DIR) -L$(TABIX_DIR)
export LDLIBS := -lbam -ltabix -lz 



STRELKA_PROGRAMS := bin/strelka2 bin/strelkaSiteSimulator
STARLING_PROGRAM := bin/starling2

STARKA_PROGRAMS := $(STRELKA_PROGRAMS) $(STARLING_PROGRAM)


export SHARED_CXX_DIR := starling_common common blt_common blt_util
export SHARED_CXX_LIB  := $(foreach dir,$(SHARED_CXX_DIR),$(CURDIR)/lib/$(dir)/lib$(dir).a)

STRELKA_OBJS := $(patsubst %.cpp,%.o, $(wildcard lib/strelka/*.cpp))
STAR_OBJS := $(patsubst %.cpp,%.o, $(wildcard lib/starling/*.cpp))


.PHONY: all build clean install test

# not really phony, but we want to force a check of whether source files are newer than lib:
.PHONY: $(SHARED_CXX_LIB)


all: build

build: $(STARKA_PROGRAMS)

$(SHARED_CXX_LIB) :
	$(MAKE) -C $(dir $@) $(notdir $@)

install: build
	cp $(STARKA_PROGRAMS) $(BIN_DIR)

test: build cppcheck.done
	$(MAKE) -C lib $@

cppcheck.done: $(STARKA_PROGRAMS)
	@echo "Running cppcheck" && \
	$(CURDIR)/run_cppcheck.bash

$(STARKA_PROGRAMS): LDLIBS += -lboost_program_options
$(STARKA_PROGRAMS): CC := $(CXX) # this fixes the linker to inlcude libc++

$(STRELKA_PROGRAMS): % : %.o $(STRELKA_OBJS) $(SHARED_CXX_LIB)
$(STARLING_PROGRAM): % : %.o $(STAR_OBJS) $(SHARED_CXX_LIB)


clean:
	$(MAKE) -C lib $@
	$(RM) $(STRELKA_OBJS) $(STAR_OBJS) $(SHARED_CXX_LIB) $(STARKA_PROGRAMS)


# extra dev tags:
#
etags: TAGS

TAGS:
	find . -type f -iname "*.[ch]" -or -iname "*.cpp" -or -iname "*.hh" | sed "s/\.\///" | etags - 


