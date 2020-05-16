// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for SIR on temporal networks by Petter Holme (2018)

// include header file
#include "main.h"

// declare g as a GLOBALS struct
GLOBALS g;

// declare n as a array of NODE structs
NODE *n;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine first localizes the first contact later than 'now' in t
// then picks the contact that can infect (a chain of Bernoulli trials)
// among the rest of the contacts. It returns the time of the infecting contact
// *t should be an array of contact times, nt should be the number of contacts, 
// and now is the current time

unsigned int next_contact (unsigned int *t, unsigned int nt, unsigned int now) {
	
    // declare some integers
    // hi = nt - 1 because index starts at 0
    unsigned int i, lo = 0, mid, hi = nt - 1;

    // no need to search further because t is sorted
    // return END which is max unsigned integer
	if (t[hi] <= now) return END;

	// the actual bisection search
	do {
        // >> is a right shift operator (here bits of '(lo + hi)' are shifted to the right by 1 place)
        // effectively, this computes floor((lo + hi) / 2)
		mid = (lo + hi) >> 1;
		if (t[mid] > now) hi = mid;
		else lo = mid;
	} while (hi > lo + 1);

    // the only case lo is correct
    // if now is smaller than first contact time
	if (now < t[lo]) hi = lo;

	// get a random contact
	i = hi + g.rnd2inx[pcg_16()];

    // if the resulting contact is too late, skip it
    // integer i is larger or equal the number of contacts
    // NONE is max. unsigned integer - 1
	if (i >= nt) return NONE;

	// return the time of the contact
	return t[i];
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine does the book keeping for an infection event

void infect () {
    
    // declare integers, me is first element of heap
	unsigned int i, you, t, me = g.heap[1];
    
    // declare integers, now is the infection time of me, duration is recovery time of me
	unsigned int now = n[me].time, duration = exptime();

    // take the newly infected off the heap
	del_root();

    // if the duration is zero, no one else can be infected
	if (duration > 0) {
        
        // set time of me to when me recovers
		n[me].time += duration;

		// go through the neighbors of the infected node
		for (i = 0; i < n[me].deg; i++) {
            
            // get neighbor i as you
			you = n[me].nb[i];
            
            // if you is S (susceptible), you can be infected
            // S is a macro defined in header file
			if (S(you)) {
                
				// find the infection time of you
				t = next_contact(n[me].t[i], n[me].nc[i], now);
                
                // bcoz the sorting of nbs, we can break
                // hmmm, so neighbors seem to be sorted so that we don't need to continue for loop
				if (t == END) break;

				// if the infection time of you is before when me gets recovered,
				// and (if it was already listed for infection) before the
				// previously listed infection event, then list it
				if ((t <= n[me].time) && (t < n[you].time)) {
					
                    // set you's infection time
                    n[you].time = t;
                    
                    // if not listed before, then extend the heap
					if (n[you].heap == NONE) {
						g.heap[++g.nheap] = you;
						n[you].heap = g.nheap;
					}
                    
                    // this works because there the only heap relationship that can be 
                    // violated is the one between you and its parent
                    // here, we upheap if necessary
					up_heap(n[you].heap);
				}
			}
		}
	}

    // to get the outbreak size
    // we add me to the array g.s that stores all node indices that have been infected
    // no matter whether or not they recovered (infected or recovered is not differentiated)
	g.s[g.ns++] = me;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine runs one SIR outbreak from a random starting node

void sir () {
    
    // declare unsigned integers
	unsigned int i, source;
	
    // initialize size of outbreak to 0
	g.ns = 0;
	
    // randomly select the source of the outbreak
	source = pcg_32_bounded(g.n);
    // randomly select the time of the infection for source
	n[source].time = pcg_32_bounded(g.dur);
    
    // set heap of source to 1
	n[source].heap = 1;
    // put source at index 1 of heap
	g.heap[g.nheap = 1] = source;

	// run the outbreak as long as heap contains elements
	while (g.nheap) infect();

	// set time of infected and recovered nodes back to NONE
	for (i = 0; i < g.ns; i++) {
        // increase ni for all nodes in g.s
        n[g.s[i]].ni++;
        // set heap and time back to NONE
        n[g.s[i]].heap = n[g.s[i]].time = NONE;
    }
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// funs the full simulation procedure

void simulate () {
    
    // declare integers i
	unsigned int i;
    
    // allocate memory to g.s (array containing the nodes that are infected/recovered
    g.s = calloc(g.n, sizeof(unsigned int));

	// initialize so that heap and time of every node are set to infinity (or NONE)
	for (i = 0; i < g.n; i++) {
        // initialize all ni to 0
        n[i].ni = 0;
        // set heap and time for every node to NONE
        n[i].heap = n[i].time = NONE;
    }
	
	// run the simulations
	for (i = 0; i < NSIM; i++) {
        
        // run sir() NSIM times
		sir();
        
        // print progress bar
        progress_bar("Simulation progress: ", i, NSIM);

	}
    
    // free memory allocated to g.s
	free(g.s);

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
