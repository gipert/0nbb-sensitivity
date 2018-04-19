# Makefile
#
# Author: Luigi Pertoldi - pertoldi@pd.infn.it
# Created: 12 Apr 2018
#
EXE      = discSensVsBI
CXX      = c++
CXXFLAGS = -std=c++1y -Wall -O3 \
           $$(root-config --cflags) \
           $$(bat-config --cflags) \
           -Iinclude -Itools/jsoncpp -Itools/progressbar
ifeq ($(shell uname -s),Darwin)
	CXXFLAGS += -Wno-unused-command-line-argument
else
	CXXFLAGS += -fopenmp -DGOPARALLEL
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

bin/$(EXE) : $(OBJECTS) obj/jsoncpp.o obj/progressbar.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(CXXLIBS)

.SECONDEXPANSION:
obj/%.o : src/%.cxx $$(wildcard include/%.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(CXXLIBS)

# third-party
obj/jsoncpp.o : $(shell find tools/jsoncpp/* -type f)
	$(CXX) $(CXXFLAGS) -c -o $@ tools/jsoncpp/jsoncpp.cpp

obj/progressbar.o : $(shell find tools/progressbar/* -type f)
	$(CXX) $(CXXFLAGS) -c -o $@ tools/progressbar/ProgressBar.cxx

.PHONY : clean
clean :
	-rm -rf bin obj
