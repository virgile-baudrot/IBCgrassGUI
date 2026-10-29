# Build the IBC-grass command-line model (static binary, as in setup.sh and Dockerfile)
#   make -C Model-files        -> Model-files/IBCgrassGUI
CXX      ?= g++
CXXFLAGS ?= -std=c++11 -O2
LDFLAGS  ?= -static -static-libgcc -static-libstdc++

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

IBCgrassGUI: $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(OBJ)

%.o: %.cpp $(wildcard *.h)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) IBCgrassGUI

.PHONY: clean
