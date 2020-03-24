# Declaring the compiler, the required flags and dependencies
CC = mpicxx
CXXFLAGS = -std=c++11 -Wall -O2 -pedantic
LDLIBS = -lboost_program_options -llapack -lblas -lscalapack-openmpi
HDRS = LidDrivenCavity.h Poisson2DSolver.h
TARGET = Solve
OBJS = LidDrivenCavitySolver.o LidDrivenCavity.o Poisson2DSolver.o

# defining the make target
default: $(TARGET)
all: $(TARGET)

# generating the object files
%.o : %.cpp $(HDRS)
	$(CC) $(CXXFLAGS) -o $@ -c $< $(LDLIBS)

# linking the object file and producing the executable
$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

.PHONY: test testSerial clean cleanImages cleanData cleanAll

# the parallel test case for the code
test:$(TARGET)
	mpiexec -np 4 ./$(TARGET) --Nx=5 --Ny=5 --Px=2 --Py=2

# the serial test case for the code
testSerial:$(TARGET)
	mpiexec -np 1 ./$(TARGET) --Nx=5 --Ny=5

# cleaning commands
clean:
	rm -f $(TARGET) *.o

cleanImages:
	rm -f Images/*

cleanData:
	rm -f Data/*

cleanAll:
	rm -f Data/* Images/* *.o $(TARGET)
