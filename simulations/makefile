CC=icc -std=c++11
#CFLAGS=-c -Wall -Wextra
SOURCES= atom.cpp modelstats.cpp lattice.cpp  main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MagStr
BIN=./bin2/

$(EXECUTABLE): $(SOURCES)
	$(CC)  $(SOURCES) -o $(BIN)$(EXECUTABLE)

#$(EXECUTABLE): $(OBJECTS)
#	$(CC) $(OBJECTS) -o $@
