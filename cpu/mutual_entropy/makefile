
CC = icc
FLAGS = -O3 -g -qopenmp
BIN = epi
OBJ = $(BIN).o


$(BIN): $(OBJ)
	$(CC) -o $@ $^ $(FLAGS)

%.o: %.cpp
	$(CC) -c -o $@ $< $(FLAGS)

clean:
	rm *.o
	rm $(BIN)
