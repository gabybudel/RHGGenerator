CXX = g++

CFLAGS  = -Wall -g -Os -std=c++11
INCLUDES = -I/usr/local/include
LDFLAGS = -lm -lgsl -lgslcblas -L/usr/local/lib
SRC = main.cpp
OBJ = main.o
TARGET = generate_rhg
	
$(TARGET): $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $< $(LDFLAGS)	
	
$(OBJ): $(SRC)
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(OBJ)