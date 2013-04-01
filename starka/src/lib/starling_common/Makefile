#
# IDS must be defined:
#
IDS := CXXFLAGS LDFLAGS LDLIBS SHARED_CXX_LIB

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(IDS), $(call assert-defined-indirect,$I))

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=%.o)

TEST_SRCS = $(wildcard test/*.cpp)
TEST_OBJS = $(TEST_SRCS:%.cpp=%.o)

LIBLABEL=lib$(notdir $(CURDIR))
LIBNAME=$(LIBLABEL).a

.PHONY: clean test


$(LIBNAME): $(OBJS)
	$(AR) -csru $@ $(OBJS)

test: test_execute

TEST_PROGRAM := $(LIBLABEL)_unittest

test_execute: $(TEST_PROGRAM)
	@echo Unit testing $(LIBNAME)
	@$(TEST_PROGRAM) --log_level=test_suite

$(TEST_PROGRAM): LDLIBS += -lboost_unit_test_framework
$(TEST_PROGRAM): $(LIBNAME) $(TEST_OBJS)
	$(CXX) $(TEST_OBJS) $(SHARED_CXX_LIB) -o $(TEST_PROGRAM) $(LDFLAGS) $(LDLIBS)

$(TEST_OBJS): CXXFLAGS += -I$(CURDIR)

clean:
	rm -f $(LIBNAME) $(OBJS) $(TEST_PROGRAM) $(TEST_OBJS)

