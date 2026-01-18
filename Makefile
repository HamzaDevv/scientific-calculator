# Zillion Calculator Makefile

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
TARGET = zillion
SRC = zillion_calculator_complete.cpp

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)

# Debug build
debug: CXXFLAGS = -std=c++17 -Wall -Wextra -g -DDEBUG
debug: clean $(TARGET)
