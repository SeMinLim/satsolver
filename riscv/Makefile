LIB = -lrt

all: $(wildcard *.cpp)
	mkdir -p obj
	riscv64-unknown-elf-g++ -o obj/main $(wildcard *.cpp) -Wall -pedantic -g -O3
