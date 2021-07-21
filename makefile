CXX = g++

CFLAGS  = -Wall -g -Os -std=c++11
INCLUDES = -I/usr/local/include
LDFLAGS = -lm -lgsl -lgslcblas -L/usr/local/lib
TARGET = generate_rhg
	
$(TARGET): main.o
	$(CXX) $(CFLAGS) -o $(TARGET) main.o $(LDFLAGS)	
	
main.o: main.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c main.cpp -o main.o

.PHONY: clean
clean:
	rm -rf main.o