# Include files
SOURCES=initLB.c visualLB.c boundary.c collision.c streaming.c computeCellValues.c main.c helper.c

# Compiler
# --------
CC=gcc

CFLAGS=-Werror -pedantic -Wall -std=c99 -O3

# Linker flags
# ------------
LDFLAGS=-lm -lrt

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=lbsim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

clean:
	rm -f $(OBJECTS) $(EXECUTABLE) *.vtk


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
