CC = gcc
FILES = main.c utilityFunctions.c
OBJECTS = main.o utilityFunctions.o

main.run: $(OBJECTS)
	$(CC) -Llib $(OBJECTS) -lgsl -lgslcblas -lm

$(OBJECTS): $(FILES)
	$(CC) -Wall -Iinclude -c $(FILES)
