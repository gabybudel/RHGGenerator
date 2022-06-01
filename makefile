CXX = g++

CFLAGS  = -Wall -g -Os -std=c++11
INCLUDES = -I/usr/local/include -I/opt/homebrew/include
LDFLAGS = -lm -lgsl -lgslcblas -L/usr/local/lib -L/opt/homebrew/lib
SRC = main.cpp
OBJ = main.o
AUX = lambert
TARGET = generate_rhg
	
$(TARGET): $(OBJ) $(AUX).o
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)	
	
$(OBJ): $(SRC) $(AUX).h
	$(CXX) $(CFLAGS) $(INCLUDES) -c $(SRC)
	
$(AUX).o: $(AUX).h

.PHONY: clean

clean:
	rm -f $(OBJ) $(AUX).o