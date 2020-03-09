CC = g++
CXXFLAGS = -std=c++11 -Wall -O2 -pedantic
LDLIBS = -lboost_program_options -llapack -lblas
HDRS = LidDrivenCavity.h
TARGET = Solve
OBJS = LidDrivenCavitySolver.o LidDrivenCavity.o

default: $(TARGET)
all: $(TARGET)

%.o : %.cpp $(HDRS)
	$(CC) $(CXXFLAGS) -o $@ -c $< $(LDLIBS)

$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

.PHONY: clean run

# An example of using the code
run: $(TARGET)
	./$(TARGET)

test: $(TARGET)
	./$(TARGET)
	./$(TARGET) --Lx=1.0 --Ly=1.0 --Nx=161 --Ny=161 --Px=0 --Py=0 --dt=0.0001 --T=1.0 --Re=100
	./$(TARGET) --Lx=1.0 --Ly=1.0 --Nx=161 --Ny=161 --Px=0 --Py=0 --dt=0.0001 --T=1.0 --Re=400
	./$(TARGET) --Lx=1.0 --Ly=1.0 --Nx=161 --Ny=161 --Px=0 --Py=0 --dt=0.0001 --T=1.0 --Re=1000
	./$(TARGET) --Lx=1.0 --Ly=1.0 --Nx=161 --Ny=161 --Px=0 --Py=0 --dt=0.0001 --T=1.0 --Re=3200
	./$(TARGET) --Lx=1.0 --Ly=2.0 --Nx=161 --Ny=161 --Px=0 --Py=0 --dt=0.0001 --T=1.0 --Re=100
	./$(TARGET) --Lx=2.0 --Ly=1.0 --Nx=161 --Ny=161 --Px=0 --Py=0 --dt=0.0001 --T=1.0 --Re=100
	
clean:
	rm -f $(TARGET) *.o
