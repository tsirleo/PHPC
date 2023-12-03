CFLAGS_P = -O3 -Wall -Wextra -std=c++11
CC_P = mpic++

run_polus: compile
	for N in 32768 262144 1048576 16776256 134217728 536870912 1073741824 17179869184 ; do \
		for p in 1 2 4 8 16 32 48 56 ; do \
			mpisubmit.pl -p $$p -w 00:30 --stdout attempt\_$$p\_$$N.out --stderr attempt\_$$p\_$$N.err ./pbsort $$N ; \
		done \
	done \

run: compile
	mpiexec ./pbsort 100000

compile: clean
	$(CC_P) $(CFLAGS_P) parallel_bsort.cpp -o pbsort

clean:
	rm -f *.o pbsort
