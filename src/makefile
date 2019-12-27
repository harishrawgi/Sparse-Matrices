a.out : main.o sparseMatrix.o utils.o
	gcc main.o sparseMatrix.o utils.o

main.o: main.c sparseMatrix.h
	gcc -c main.c

sparseMatrix.o : sparseMatrix.c sparseMatrix.h utils.h
	gcc -c sparseMatrix.c

utils.o : utils.c utils.h
	gcc -c utils.c
