OBJ = main.o CellModel.o CellModelSetup.o

CC = g++
CFLAGS = -g # -Wall
# INC = -I/usr/local/boost_1_47_0
# LIBS = -lpthread
LIBS += -lncurses

all: celldiff

CellModel.o: CellModel.cpp
	$(CC) $(CFLAGS) $(INC) -c -o CellModel.o CellModel.cpp 
	
CellModelSetup.o: CellModelSetup.cpp
	$(CC) $(CFLAGS) $(INC) -c -o CellModelSetup.o CellModelSetup.cpp 
	
main.o: main.cpp
	$(CC) $(CFLAGS) $(INC) -c -o main.o main.cpp

celldiff: $(OBJ)
	$(CC) $(CFLAGS) $(INC) $(OBJ) -o celldiff $(LIBS)

clean:
	rm *.o
	rm celldiff
	
.PHONY: all clean
