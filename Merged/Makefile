CXX := g++
CXXFLAGS := -std=c++17 $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
BATLIBS := -L/usr/local/lib -lBAT
SRC := runfit.cpp bkgModel.cpp sigModel.cpp
OBJ := $(SRC:.cpp=.o)
TARGET := runfit

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS) $(BATLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -Wall -Wextra -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
