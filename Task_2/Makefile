CFLAGS_P = -O3 -Wall -Wextra -std=c++11 -fopenmp
CC_P = mpic++

run_polus: compile
	for N in 1024 8192 65536 524288 4194304 33554432 67108864 134217728 ; do \
 		for p in 1 2 4 8 16 32 48 64 ; do \
 			mpisubmit.pl -p $$p -t 8 -w 00:30 --stdout attempt\_$$p\_$$N.out --stderr attempt\_$$p\_$$N.err "./pbsort --size=$$N" ; \
 		done \
	done \

run: compile
	mpiexec ./pbsort --size=100000

compile: clean
	$(CC_P) $(CFLAGS_P) parallel_bsort.cpp -o pbsort

clean:
	rm -f *.o pbsor
