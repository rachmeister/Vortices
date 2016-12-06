CC=gcc
CFLAGS=-Wall -Werror
DEPS = def.h methods.h
LDFLAGS = -lm
OBJS = main.o util.o integrators.o test_cases.o 

all: vortices

vortices: $(OBJS)
	$(CC) $(LDFLAGS) $(CFLAGS) $(OBJS) -o vortices

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm -f *~ *.o vortices

