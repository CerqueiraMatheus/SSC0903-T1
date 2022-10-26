seq:
	gcc studentsseq.c -g -Wall -Werror -lm -fopenmp -Ofast

par:
	gcc studentspar.c -g -Wall -Werror -lm -fopenmp -Ofast

run:
	./a.out