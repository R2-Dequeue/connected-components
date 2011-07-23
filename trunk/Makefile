NAME = ConnectedComponents

SOURCES = algebraic.cpp cad.cpp main.cpp polynomialbase.cpp polynomialq.cpp polynomialqq.cpp
OBJECTS = $(SOURCES:.cpp=.o)

CC = mingw32-g++
RM = del

BOOSTDIR = C:\dev\boost_1_45_0
INCLUDEDIR = C:\MinGW\msys\1.0\local\include
LIBDIR = C:\MinGW\msys\1.0\local\lib

DEBUG = -DNDEBUG -s
LIBS = -lginac -lcln -lgmp

CFLAGS = -Wall -fexceptions -O2 $(DEBUG) -I$(BOOSTDIR) -I$(INCLUDEDIR)
LFLAGS = $(DEBUG) -L$(LIBDIR)

# Note: Do NOT put LIBs up here; leave it below.

###############################################################################

$(NAME).exe: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o $@ $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o $(NAME).exe

#	$@ name of the target
#	$^ name of all prerequisites with duplicates removed
#	$< name of the first prerequisite
