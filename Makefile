CC = g++
CXXFLAGS = -Wall -Werror -O3 $(shell pkg-config --cflags jellyfish-2.0) -std=c++11
# this is very unportable
LDFLAGS := $(LDFLAGS) $(shell pkg-config --libs-only-L jellyfish-2.0) -L$(BOOST_ROOT)/lib 
LDLIBS = -lboost_timer -lboost_chrono -lboost_system -lboost_program_options  $(shell pkg-config --libs-only-l jellyfish-2.0)

all: bin/jellyfish-matrix bin/jellyfish-matrix-histo bin/jellyfish-matrix-read
bin/jellyfish-matrix: jellyfish-matrix
	cp $< $@
bin/jellyfish-matrix-histo: jellyfish-matrix-histo
	cp $< $@
bin/jellyfish-matrix-read: jellyfish-matrix-read
	cp $< $@

jellyfish-matrix: jellyfish-matrix.cc
jellyfish-matrix-histo: jellyfish-matrix-histo.cc
jellyfish-matrix-read: jellyfish-matrix-read.cc
clean:
	rm -f bin/* jellyfish-matrix jellyfish-matrix-histo jellyfish-matrix-read
