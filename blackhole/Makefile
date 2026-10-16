CC ?= gcc
CFLAGS ?= -O3 -fopenmp -Wall -Wextra

blackhole: blackhole.c
	$(CC) $(CFLAGS) -o $@ $< -lm

clean:
	rm -f blackhole *.ppm
