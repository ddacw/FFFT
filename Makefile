CXX = g++
CFLAGS = -pthread -ltbb -std=c++17 -Wall

SOURCES = main.cpp pfft.h refft.h main.h

main: main.cpp
	$(CXX) $(CFLAGS) -o main main.cpp

clean:
	rm -f main
