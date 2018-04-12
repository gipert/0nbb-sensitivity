# Makefile
#
# Author: Luigi Pertoldi - pertoldi@pd.infn.it
# Created: 12 Apr 2018
#
EXE      = discSensVsBI
CXX      = c++
CXXFLAGS = -std=c++0y -Wall \
           $$(root-config --cflags) \
           $$(bat-config --cflags) \
           -Iinclude
ifeq ($(shell uname -s),Darwin)
	CXXFLAGS += -Wno-unused-command-line-argument
endif
CXXLIBS  = $$(root-config --libs) \
           $$(bat-config --libs)
DIRS     = bin obj

.PHONY : all
all : $(DIRS) bin/$(EXE)

$(DIRS) :
	mkdir -p $@

SOURCES = $(wildcard src/*.cxx)
OBJECTS = $(patsubst src/%.cxx, obj/%.o, $(SOURCES))

bin/$(EXE) : $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(CXXLIBS)

.SECONDEXPANSION:
obj/%.o : src/%.cxx $$(wildcard include/%.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(CXXLIBS)

.PHONY : clean
clean :
	-rm -rf bin obj
