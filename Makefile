.PHONY: clean
eigen: main.o eigen.o schur.o hessenberg.o complex_matrix.o
	gcc -o eigen -Wall main.o eigen.o schur.o hessenberg.o complex_matrix.o -lm
main.o: main.c
	gcc -c -Wall main.c
eigen.o: eigen.c
	gcc -c -Wall eigen.c
schur.o: schur.c
	gcc -c -Wall schur.c
hessenberg.o: hessenberg.c
	gcc -c -Wall hessenberg.c
complex_matrix.o: complex_matrix.c
	gcc -c -Wall complex_matrix.c
clean:
	rm *.o
