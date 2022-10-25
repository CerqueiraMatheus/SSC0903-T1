seq:
	gcc studentsseq.c -g -Wall -Werror -lm -o seq.out

par:
	gcc studentspar.c -fopenmp -lm -g -Wall -Werror -o par.out