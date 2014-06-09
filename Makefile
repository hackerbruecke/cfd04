# Include files
SOURCES=initLB.c visualLB.c boundary.c collision.c streaming.c computeCellValues.c main.c helper.c

# Compiler
# --------
CC=mpicc

CFLAGS=-Werror -pedantic -Wall -std=c99 -O0

# Linker flags
# ------------
LDFLAGS=-lm -lrt

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=lbsim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

clean:
	rm -f $(OBJECTS) $(EXECUTABLE) *.vtk *.pvts


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
