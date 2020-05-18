SRC = .
CFLAGS = -W -Wall -DTIME -Ofast -march=native
LDFLAGS = 
CC = gcc

OBJ1 = o/main.o o/tsir.o o/misc.o o/heap.o o/quick.o o/pcg_rnd.o o/greedy.o

all : main

main: $(OBJ1)
	$(CC) $(LDFLAGS) -o $@ $^ -lm

o/main.o : $(SRC)/main.c $(SRC)/main.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/main.c -o $@

o/tsir.o : $(SRC)/tsir.c $(SRC)/main.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/tsir.c -o $@

o/misc.o : $(SRC)/misc.c $(SRC)/main.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/misc.c -o $@

o/heap.o : $(SRC)/heap.c $(SRC)/main.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/heap.c -o $@

o/quick.o : $(SRC)/quick.c $(SRC)/main.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/quick.c -o $@

o/pcg_rnd.o : $(SRC)/pcg_rnd.c $(SRC)/main.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/pcg_rnd.c -o $@
    
o/greedy.o : $(SRC)/greedy.c $(SRC)/main.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/greedy.c -o $@
