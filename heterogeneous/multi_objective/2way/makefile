
CC = icc
FLAGS = -I/opt/intel/inteloneapi/compiler/latest/linux/lib/oclfpga/host/include/ -lOpenCL -O3 -g -std=c++11 -qopenmp
BIN = epi_2way
OBJ = $(BIN).o


$(BIN): $(OBJ)
	$(CC) -o $@ $^ $(FLAGS)
	
%.o: %.cpp
	$(CC) -c -o $@ $< $(FLAGS)

clean:
	rm *.o
	rm $(BIN)
