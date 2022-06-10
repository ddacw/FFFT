CXX = g++
CFLAGS = -Iinclude -pthread -ltbb -std=c++17 -Wall

SOURCES = main.cpp test.cpp src/fft.cpp src/util.cpp
OBJECTS = util.o fft.o 

test: test.cpp $(OBJECTS)
	$(CXX) $(CFLAGS) -o test $(OBJECTS) test.cpp

main: main.cpp $(OBJECTS)
	$(CXX) $(CFLAGS) -o main $(OBJECTS) main.cpp

fft.o: src/fft.cpp
	$(CXX) -c $(CFLAGS) -o fft.o src/fft.cpp

util.o: src/util.cpp
	$(CXX) -c $(CFLAGS) -o util.o src/util.cpp

clean:
	rm -f *.o
	rm -f main test
	rm -f fft_*.txt
