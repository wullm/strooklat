#Compiler options
GCC = gcc

CFLAGS = -Wall -fopenmp -march=native -O4 -Wshadow -fPIC

all:
	$(GCC) strooklat.c -c -o strooklat.o $(CFLAGS)
	$(GCC) test.c -o test strooklat.o $(CFLAGS) -lm
	$(GCC) test_inline.c -o test_inline $(CFLAGS) -lm

format:
	clang-format-10 -style="{BasedOnStyle: LLVM, IndentWidth: 4, AlignConsecutiveMacros: true, IndentPPDirectives: AfterHash}" -i *.c *.h
	astyle -i *.h *.c

check:
	cd tests && make
