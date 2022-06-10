CXX = g++
CFLAGS = -pthread -ltbb -std=c++17 -Wall

SOURCES = main.cpp main.h fft.h gfft.h utils.h
OBJECTS = utils.o fft.o main.o 

main: $(OBJECTS)
	$(CXX) $(CFLAGS) -o main $(OBJECTS)

utils.o: utils.cpp utils.h
	$(CXX) -c $(CFLAGS) -o utils.o utils.cpp

fft.o: fft.cpp fft.h utils.h
	$(CXX) -c $(CFLAGS) -o fft.o fft.cpp

main.o: main.cpp main.h
	$(CXX) -c $(CFLAGS) -o main.o main.cpp

clean:
	rm -f *.o
	rm -f main
