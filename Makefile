all: wifi wall

wall: farwall-greedy.cpp
	g++ farwall-greedy.cpp -o wall -std=gnu++11 -O2

wifi: main.cpp
	g++ main.cpp -o wifi -std=gnu++11 -O2
