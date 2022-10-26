seq:
	gcc studentsseq.c -g -Wall -Werror -lm -fopenmp -O3

par:
	gcc studentspar.c -g -Wall -Werror -lm -fopenmp -O3

run:
	./a.out