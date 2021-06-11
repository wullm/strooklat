#Compiler options
GCC = gcc

LIBRARIES = -lgsl -lgslcblas
CFLAGS = -Wall -fopenmp -march=native -O4 -Wshadow -fPIC

all:
	# $(GCC) strooklat.c -c -o strooklat.o $(CFLAGS)
	# $(GCC) test.c -o test strooklat.o $(CFLAGS) -lm $(LIBRARIES)
	$(GCC) test.c -o test $(CFLAGS) -lm $(LIBRARIES)
	$(GCC) test_nd.c -o test_nd $(CFLAGS) -lm $(LIBRARIES)

format:
	clang-format-10 -style="{BasedOnStyle: LLVM, IndentWidth: 4, AlignConsecutiveMacros: true, IndentPPDirectives: AfterHash}" -i *.c *.h
	astyle -i *.h *.c

check:
	cd tests && make
