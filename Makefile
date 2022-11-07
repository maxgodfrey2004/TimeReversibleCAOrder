.PHONY: data run all

REXE := Rscript.exe
CC := g++

DATA_FILE_NAME := aadata.txt

all: automata.out 

run: automata.out
	./automata.out 10 10

automata.out: automata.cpp
	$(CC) automata.cpp -o automata.out
