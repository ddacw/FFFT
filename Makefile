CXX = g++
CFLAGS = -Iinclude -pthread -ltbb -std=c++17 -Wall

SOURCES = main.cpp test.cpp src/mfft.cpp src/fft.cpp src/util.cpp
OBJECTS = mfft.o fft.o util.o

.PHONY : all

all : main test

main: main.cpp $(OBJECTS)
	$(CXX) $(CFLAGS) -o main $(OBJECTS) main.cpp

test: test.cpp $(OBJECTS)
	$(CXX) $(CFLAGS) -o test $(OBJECTS) test.cpp

mfft.o: src/mfft.cpp
	$(CXX) -c $(CFLAGS) -o mfft.o src/mfft.cpp

fft.o: src/fft.cpp
	$(CXX) -c $(CFLAGS) -o fft.o src/fft.cpp

util.o: src/util.cpp
	$(CXX) -c $(CFLAGS) -o util.o src/util.cpp

clean:
	rm -f *.o
	rm -f main test
	rm -f fft_*.txt
