CC = mpicxx
CFLAGS = -g
LDFLAGS = -L$(HOME)/sw/lib
INCLUDE = -I$(HOME)/sw/include
LIB = -lntl -lmpc -lmpfr -lgmp -lm
OUTPUT_A2M = -o $(EXECUTABLE_A2M)
EXECUTABLE_A2M = A2Mult
OUTPUT_A2A = -o $(EXECUTABLE_A2A)
EXECUTABLE_A2A = A2Add
OUTPUT_B2M = -o $(EXECUTABLE_B2M)
EXECUTABLE_B2M = B2Mult
OUTPUT_G2M = -o $(EXECUTABLE_G2M)
EXECUTABLE_G2M = G2Mult
OUTPUT_G2A = -o $(EXECUTABLE_G2A)
EXECUTABLE_G2A = G2Add
OUTPUT_HasseM = -o $(EXECUTABLE_HasseM)
EXECUTABLE_HasseM = HessianMult
OUTPUT_HasseA = -o $(EXECUTABLE_HasseA)
EXECUTABLE_HasseA = HessianAdd
OUTPUT_K3 = -o $(EXECUTABLE_K3)
EXECUTABLE_K3 = Kloosterman3
OUTPUT_A3M = -o $(EXECUTABLE_A3M)
EXECUTABLE_A3M = A3Mult
OUTPUT_B3M = -o $(EXECUTABLE_B3M)
EXECUTABLE_B3M = B3Mult
RUNTIME_LIB = -Wl,-rpath=$(HOME)/sw/lib
OPENMP = -openmp


default:
	$(CC) $(INCLUDE) A2Mult.cxx $(OUTPUT_A2M) $(LDFLAGS) $(LIB) $(RUNTIME_LIB) $(OPENMP)
	$(CC) $(INCLUDE) A2Add.cxx $(OUTPUT_A2A) $(LDFLAGS) $(LIB) $(RUNTIME_LIB) $(OPENMP)
	$(CC) $(INCLUDE) B2Mult.cxx $(OUTPUT_B2M) $(LDFLAGS) $(LIB) $(RUNTIME_LIB) $(OPENMP)
	$(CC) $(INCLUDE) G2Mult.cxx $(OUTPUT_G2M) $(LDFLAGS) $(LIB) $(RUNTIME_LIB) $(OPENMP)
	$(CC) $(INCLUDE) G2Add.cxx $(OUTPUT_G2A) $(LDFLAGS) $(LIB) $(RUNTIME_LIB) $(OPENMP)
	$(CC) $(INCLUDE) HessianMult.cxx $(OUTPUT_HasseM) $(LDFLAGS) $(LIB) $(RUNTIME_LIB) $(OPENMP)
	$(CC) $(INCLUDE) HessianAdd.cxx $(OUTPUT_HasseA) $(LDFLAGS) $(LIB) $(RUNTIME_LIB) $(OPENMP)
	$(CC) $(INCLUDE) Kloosterman3.cxx $(OUTPUT_K3) $(LDFLAGS) $(LIB) $(RUNTIME_LIB) $(OPENMP)
	$(CC) $(INCLUDE) A3Mult.cxx $(OUTPUT_A3M) $(LDFLAGS) $(LIB) $(RUNTIME_LIB) $(OPENMP)
	$(CC) $(INCLUDE) B3Mult.cxx $(OUTPUT_B3M) $(LDFLAGS) $(LIB) $(RUNTIME_LIB) $(OPENMP)

clean:
	rm $(EXECUTABLE_A2M)
	rm $(EXECUTABLE_A2A)
	rm $(EXECUTABLE_B2M)
	rm $(EXECUTABLE_G2M)
	rm $(EXECUTABLE_G2A)
	rm $(EXECUTABLE_HasseM)
	rm $(EXECUTABLE_HasseA)
	rm $(EXECUTABLE_K3)
	rm $(EXECUTABLE_A3M)
	rm $(EXECUTABLE_B3M)