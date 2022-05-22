CXX = g++
CFLAGS = -pthread -std=c++17 -Wall

SOURCES = main.cpp 

main: main.cpp
	$(CXX) $(CFLAGS) -o main main.cpp

clean:
	rm -f main
