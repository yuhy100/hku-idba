SOURCES:=$(wildcard *.cpp)
OBJECTS:=$(SOURCES:%.cpp=%.o)
TARGETS:=$(SOURCES:%.cpp=%)
LIBS:=../libs/libs.a -lpthread

CC:=g++
CXXFLAGS=-Wall -I../libs -O2
LDFLAGS=$(LIBS)

all:$(TARGETS)
	cp $(TARGETS) ../../bin

$(TARGETS):%:%.o $(LIBS)
	$(CC) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGETS) $(OBJECTS) *~

