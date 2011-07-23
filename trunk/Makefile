# Comment out the DEBUG macro to turn off debugging symbols.

NAME = ConnectedComponents
SOURCES = algebraic.cpp cad.cpp main.cpp polynomialbase.cpp polynomialq.cpp polynomialqq.cpp
OBJECTS = $(SOURCES:.cpp=.o)
CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
# ginac cln gmp

#CFLAGS  = -I/home/scale/g++Projects/gLib/
#LDFLAGS = -lfltk

$(NAME).exe: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $<

clean:
	rm *.o $(NAME).exe

#	$@ name of the target
#	$^ name of all prerequisites with duplicates removed
#	$< name of the first prerequisite

#algebraic.o: algebraic.cpp
#cad.o: cad.cpp
#main.o: main.cpp
#polynomialbase.o: polynomialbase.cpp
#polynomialq.o: polynomialq.cpp
#polynomialqq.o: polynomialqq.cpp
