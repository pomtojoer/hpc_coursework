CC = mpicxx
CXXFLAGS = -std=c++11 -Wall -O2 -pedantic
LDLIBS = -lboost_program_options -llapack -lblas -lscalapack-openmpi
HDRS = LidDrivenCavity.h Poisson2DSolver.h
TARGET = Solve
OBJS = LidDrivenCavitySolver.o LidDrivenCavity.o Poisson2DSolver.o

default: $(TARGET)
all: $(TARGET)

%.o : %.cpp $(HDRS)
	$(CC) $(CXXFLAGS) -o $@ -c $< $(LDLIBS)

$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

.PHONY: test testSerial clean cleanImages cleanData cleanAll

# An example of using the code
test:$(TARGET)
	mpiexec -np 4 ./$(TARGET) --Nx=5 --Ny=5 --Px=2 --Py=2

testSerial:$(TARGET)
	mpiexec -np 1 ./$(TARGET) --Nx=5 --Ny=5

clean:
	rm -f $(TARGET) *.o

cleanImages:
	rm -f Images/*

cleanData:
	rm -f Data/*

cleanAll:
	rm -f Data/* Images/* *.o $(TARGET)
