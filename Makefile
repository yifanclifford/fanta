MATLABROOT = D:/software/matlab
OBJ= main.o fanta.o graph.o ismin.o misc.o utils.o dfs.o
CC = g++ -std=c++0x
LIB= -L$(MATLABROOT)/bin/win64 -lmex -lmat -lmx
INCLUDE = -I$(MATLABROOT)/extern/include

fanta:OPTFLAGS = -m64 -O3 -s -Wno-deprecated
fanta:$(OBJ)
	$(CC) $(OPTFLAGS) -o $@ $^

debug:DEBFLAGS = -m64 -g3 -O2 -Wno-deprecated
debug:$(OBJ)
	$(CC) $(DEBFLAGS) -o $@ $^

$(OBJ):%.o:%.cpp fanta.h utils.h
	${CC} ${OPTFLAGS} -c $<

clean:
	rm -f *.o

.PHONY: debug clean
