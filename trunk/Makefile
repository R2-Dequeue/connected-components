NAME = ConnectedComponents

CORE_SOURCES = algebraic.cpp cad.cpp polynomialbase.cpp polynomialq.cpp polynomialqq.cpp

SOURCES = $(CORE_SOURCES) main.cpp
OBJECTS = $(SOURCES:.cpp=.o)
CHECK_SOURCES = $(CORE_SOURCES) check/check.cpp
CHECK_OBJECTS = $(CHECK_SOURCES:.cpp=.o)

CC = mingw32-g++
RM = del

ROOT = C:/MinGW/msys/1.0
BOOSTDIR = C:/dev/boost_1_45_0
INCLUDEDIR = $(ROOT)/local/include
LIBDIR = $(ROOT)/local/lib

DEBUG = -DNDEBUG -s
LIBS = -lginac -lcln -lgmp

CFLAGS = -Wall -fexceptions -O2 $(DEBUG) -I$(BOOSTDIR) -I$(INCLUDEDIR)
LFLAGS = $(DEBUG) -L$(LIBDIR)

# Note: Do NOT put LIBs up here; leave below.

###############################################################################

$(NAME).exe: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o $@ $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) *.o $(NAME).exe check.exe check\*.o

check: check.exe
	check.exe

check.exe: $(CHECK_OBJECTS)
	$(CC) $(LFLAGS) $(CHECK_OBJECTS) -o $@ $(LIBS)

html:
	doxygen Doxyfile

#	$@ name of the target
#	$^ name of all prerequisites with duplicates removed
#	$< name of the first prerequisite
