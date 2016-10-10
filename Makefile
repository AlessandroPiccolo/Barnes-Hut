##
# Makefile for galsim, galaxy simulation
##

CC         =  gcc
CFLAGS 	   =  -Wall
INCLUDES   =  -I/opt/X11/include
LDFLAGS    =  -L/opt/X11/lib -lX11 -lm

galsim: galsim.o graphics.o file_operations.o
	$(CC) $(CFLAGS) galsim.o graphics.o file_operations.o -o galsim $(LDFLAGS)

galsim.o: galsim.c graphics.h file_operations.h
	$(CC) -c galsim.c $(INCLUDES)

graphics.o: graphics.h graphics.c
	$(CC) -c graphics.c $(INCLUDES)

file_operations.o: file_operations.h file_operations.c
	$(CC) -c file_operations.c $(INCLUDES)

brutegalsim: brutegalsim.o graphics.o file_operations.o
	$(CC) $(CFLAGS) brutegalsim.o graphics.o file_operations.o -o brutegalsim $(LDFLAGS)

brutegalsim.o: brutegalsim.c graphics.h file_operations.h
	$(CC) -c brutegalsim.c $(INCLUDES)

compare: compare_gal_files.o file_operations.o
	$(CC) $(CFLAGS) compare_gal_files.o file_operations.o -o compare -lm

compare_gal_files.o: compare_gal_files.c file_operations.h
	$(CC) -c compare_gal_files.c

clean:
	rm -f ./galsim *.o 
