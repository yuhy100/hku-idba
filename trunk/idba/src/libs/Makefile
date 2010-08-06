SOURCES:=$(wildcard *.cpp)
OBJECTS:=$(SOURCES:%.cpp=%.o)
TARGETS:=libs.a

CC:=g++
CXXFLAGS=-Wall -O2
LDFLAGS=

all:$(TARGETS)

$(OBJECTS):%.o:%.cpp %.h globals.h
	$(CC) -c -o $@ $< $(CXXFLAGS)

$(TARGETS):%:$(OBJECTS)
	ar or $@ $(OBJECTS)

clean:
	rm -f $(TARGETS) $(OBJECTS) *~
