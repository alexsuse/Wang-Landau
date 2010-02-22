CC = g++
OBJECTS = moralWLwindows.o 
INCLUDE = -L/usr/local/lib -I/usr/local/include

all: moralWLwindows

moralWLwindows : moralWLwindows.cpp
	$(CC)  -o moralWLwindows moralWLwindows.cpp $(INCLUDE) -lgsl -lgslcblas -lm -ligraph


clean :  $(OBJECTS)
	rm $(OBJECTS)
