CC = mpicxx
CFLAGS = -g
LDFLAGS = -L$(HOME)/sw/lib
INCLUDE = -I$(HOME)/sw/include
LIB = -lntl -lmpc -lmpfr -lgmp -lm
EXECUTABLE = -o char
RUNTIME_LIB = -Wl,-rpath=$(HOME)/sw/lib
OPENMP = -openmp


default:
	$(CC) $(INCLUDE) char.c++ $(EXECUTABLE) $(LDFLAGS) $(LIB) $(RUNTIME_LIB) 
