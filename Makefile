all: TP1

TP1:  TP1.o nrio.o nralloc.o nrarith.o
	gcc -Wall -g TP1.o nrio.o nralloc.o nrarith.o -o TP1 -lm

TP1.o:  TP1.c
	gcc -Wall -g -c TP1.c -o TP1.o

nrio.o: nrio.h nrio.c
	gcc -Wall -g -c  nrio.c -o nrio.o

nralloc.o: nralloc.h nralloc.c
	gcc -Wall -g -c  nralloc.c -o nralloc.o

nrarith.o: nrarith.h nrarith.c
	gcc -Wall -g -c  nrarith.c -o nrarith.o


clean :
	rm -f TP1 *.o
