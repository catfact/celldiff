OBJ = main.o CellModel.o

CC = g++
CFLAGS = -g
BOOST = /usr/local/boost_1_47_0
# LIBS = -lpthread
LIBS = -lboost_program_options
LPATH = -L$(BOOST)/stage/lib

all: celldiff

CellModel.o: CellModel.cpp
	$(CC) $(CFLAGS) -I$(BOOST) -c -o CellModel.o CellModel.cpp
	
main.o: main.cpp
	$(CC) $(CFLAGS) -I$(BOOST) -c -o main.o main.cpp

celldiff: $(OBJ)
	$(CC) $(CFLAGS) $(INC) $(OBJ) -o celldiff $(LIBS) $(LPATH)

clean:
	rm *.o
	rm celldiff
	
.PHONY: all clean