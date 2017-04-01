all: wifi wall

wall: farwall.cpp
	g++ farwall.cpp -o wall -std=gnu++11 -O2

wifi: main.cpp
	g++ main.cpp -o wifi -std=gnu++11 -O2
