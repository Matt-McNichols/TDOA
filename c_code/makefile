CC = gcc
FILES = main.c utilityFunctions.c
OBJECTS = main.o utilityFunctions.o

main.run: $(OBJECTS)
	$(CC) -L/home/matt/mattBin/lib $(OBJECTS) -lgsl -lgslcblas -lm

$(OBJECTS): $(FILES)
	$(CC) -Wall -I/home/matt/mattBin/include -c $(FILES)
