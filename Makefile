CC = mpicxx
CFLAGS = -g
LDFLAGS = -L$(HOME)/sw/lib
INCLUDE = -I$(HOME)/sw/include
LIB = -lntl -lmpc -lmpfr -lgmp -lm
OUTPUTA2M = -o $(EXECUTABLEA2M)
EXECUTABLEA2M = A2Mult
OUTPUTA2A = -o $(EXECUTABLEA2A)
EXECUTABLEA2A = A2Add
RUNTIME_LIB = -Wl,-rpath=$(HOME)/sw/lib
OPENMP = -openmp


default:
	$(CC) $(INCLUDE) A2Mult.cxx $(OUTPUTA2M) $(LDFLAGS) $(LIB) $(RUNTIME_LIB) $(OPENMP)
	$(CC) $(INCLUDE) A2Add.cxx $(OUTPUTA2A) $(LDFLAGS) $(LIB) $(RUNTIME_LIB) $(OPENMP)

clean:
	rm $(EXECUTABLEA2M)
	rm $(EXECUTABLEA2A)