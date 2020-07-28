// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for SIR on temporal networks by Petter Holme (2018)

// include header file
#include "main.h"

// declare g as a GLOBALS struct
extern GLOBALS g;

// declare n as a array of NODE structs
extern NODE *n;

// this is an array that helps managing the size of the 'inf' arrays
unsigned int *alloc;

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
    
    // store infection run, diff. between infection time and outbreak start time, and current size of outbreak for node me
    n[me].inf[n[me].ni] = g.sim_id;
    n[me].dtime[n[me].ni] = g.t_end - now;
    n[me].dsize[n[me].ni] = g.ns + 1;
    
    // increment ni by 1
    n[me].ni++;

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
				if ((t <= n[me].time) && (t < n[you].time) && (t < g.t_end)) {
					
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
	unsigned int source;
	
    // initialize size of outbreak to 0
	g.ns = 0;
	
    // randomly select the source of the outbreak
	source = pcg_32_bounded(g.n);
    // randomly select the time of the infection for source
	n[source].time = pcg_32_bounded_ul(g.t_start, g.t_end);
    //n[source].time = rnd_bounded(g.t_start, g.t_end);
    
    // set heap of source to 1
	n[source].heap = 1;
    // put source at index 1 of heap
	g.heap[g.nheap = 1] = source;

	// run the outbreak as long as heap contains elements
	while (g.nheap) infect();
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// runs the full simulation procedure

void simulate (unsigned int t_start, unsigned int t_end) {
    
    // declare integers i
	unsigned int i, j;
    
    // declare double for avg. outbreak size
    double s = 0.0;
    
    // set start and end time in GLOBALS
    g.t_start = t_start;
    g.t_end = t_end;
    
    // allocate memory to g.outbreak_sizes
    g.outbreak_sizes = calloc(NSIM, sizeof(unsigned int));
    
    // allocate memory to g.s (array containing the nodes that are infected/recovered
    g.s = calloc(g.n, sizeof(unsigned int));
    
    // allocate memory to alloc
    alloc = calloc(g.n, sizeof(unsigned int));

	// initialize different things
	for (i = 0; i < g.n; i++) {
        // initialize all ni to 0
        n[i].ni = 0;
        // initialize alloc to 1000 for all nodes
        alloc[i] = 1000;
        // allocate memory for arrays that store simulation runs that infect node i
        n[i].inf = calloc(1000, sizeof(unsigned int));
        n[i].dtime = calloc(1000, sizeof(unsigned int));
        n[i].dsize = calloc(1000, sizeof(unsigned int));
        // set heap and time for every node to NONE
        n[i].heap = n[i].time = NONE;
    }
    
    // initialize i to 0
    i = 0;
    
    // simulate until we have NSIM scenarios
    while (i < NSIM) {
        
        // store current simulation id in globals
        g.sim_id = i;
        
        // run sir() NSIM times
		sir();
        
        // if outbreak is smaller than MIN_OUTSIZE:
        if (g.ns < MIN_OUTSIZE) {
            
            // loop over all infected nodes
            for (j = 0; j < g.ns; j++) {
                
                // set outbreak counter back to previous value
                n[g.s[j]].ni--;
                
                // set heap and time back to NONE
                n[g.s[j]].heap = n[g.s[j]].time = NONE;
            
            }
            
            // jump to next iteration of the loop
            continue;
        }
        
        // set time of infected and recovered nodes back to NONE
        for (j = 0; j < g.ns; j++) {
            
            // check if we need to allocate more memory to 'inf', 'dtime', and 'dsize'
            if (alloc[g.s[j]] == n[g.s[j]].ni) {
                // add 1000 to alloc
                alloc[g.s[j]] += 1000;
                // reallocate memory of 'inf', 'dtime', and 'dsize'
                n[g.s[j]].inf = realloc(n[g.s[j]].inf, alloc[g.s[j]] * sizeof(unsigned int));
                n[g.s[j]].dtime = realloc(n[g.s[j]].dtime, alloc[g.s[j]] * sizeof(unsigned int));
                n[g.s[j]].dsize = realloc(n[g.s[j]].dsize, alloc[g.s[j]] * sizeof(unsigned int));
            }
            
            // compute penalty reduction for size of outbreak
            n[g.s[j]].dsize[n[g.s[j]].ni - 1] = g.ns - n[g.s[j]].dsize[n[g.s[j]].ni - 1];
            
            // set heap and time back to NONE
            n[g.s[j]].heap = n[g.s[j]].time = NONE;
        }
        
        // store outbreak size
        g.outbreak_sizes[i] = g.ns;
    
        // print progress bar
        progress_bar("Simulation progress: ", i, NSIM);
        
        // increment i
        i++;
	
    }
    
    // since not all simulation runs infect every node,
    // we reallocate memory correctly
    for (i = 0; i < g.n; i++) {
        n[i].inf = realloc(n[i].inf, n[i].ni * sizeof(unsigned int));
        n[i].dtime = realloc(n[i].dtime, n[i].ni * sizeof(unsigned int));
        n[i].dsize = realloc(n[i].dsize, n[i].ni * sizeof(unsigned int));
    }
    
    // sum up all outbreak sizes
    for (i = 0; i < NSIM; i++) s += (double) g.outbreak_sizes[i];
    
    // print average outbreak size to command line
    printf("\nAverage outbreak size in training phase is %f\n", s /= NSIM);
    
    // compute median outbreak size
    // compute_median(g.outbreak_sizes, NSIM);
    
    // compute min and max outbreak size
    compute_min(g.outbreak_sizes, NSIM);
    compute_max(g.outbreak_sizes, NSIM);
    
    // free memory allocated to g.s and alloc
	free(g.s); free(g.outbreak_sizes); free(alloc);

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// create evaluation set of scenarios

void simulate_eval (unsigned int neval, unsigned int t_start, unsigned int t_end) {
    
    // declare integers i
	unsigned int i, j;
    
    // declare double for avg. outbreak size
    double s = 0.0;
    
    // set start and end time in GLOBALS
    g.t_start = t_start;
    g.t_end = t_end;
    
    // allocate memory to g.outbreak_sizes
    g.outbreak_sizes = calloc(neval, sizeof(unsigned int));
    
    // allocate memory to g.s (array containing the nodes that are infected/recovered
    g.s = calloc(g.n, sizeof(unsigned int));
    
    // allocate memory to alloc
    alloc = calloc(g.n, sizeof(unsigned int));

	// initialize different things
	for (i = 0; i < g.n; i++) {
        // initialize all ni to 0
        n[i].ni = 0;
        // initialize alloc to 1000 for all nodes
        alloc[i] = 1000;
        // allocate memory for arrays that store simulation runs that infect node i
        n[i].inf = calloc(1000, sizeof(unsigned int));
        n[i].dtime = calloc(1000, sizeof(unsigned int));
        n[i].dsize = calloc(1000, sizeof(unsigned int));
        // set heap and time for every node to NONE
        n[i].heap = n[i].time = NONE;
    }
    
    // initialize i to 0
    i = 0;
    
    // simulate until we have 'neval' scenarios
    while (i < neval) {
        
        // store current simulation id in globals
        g.sim_id = i;
        
        // run sir() neval times
		sir();
        
        // if outbreak is smaller than MIN_OUTSIZE:
        if (g.ns < MIN_OUTSIZE) {
            
            // loop over all infected nodes
            for (j = 0; j < g.ns; j++) {
                
                // set outbreak counter back to previous value
                n[g.s[j]].ni--;
                
                // set heap and time back to NONE
                n[g.s[j]].heap = n[g.s[j]].time = NONE;
            
            }
            
            // jump to next iteration of the loop
            continue;
        }
        
        // set time of infected and recovered nodes back to NONE
        for (j = 0; j < g.ns; j++) {
            
            // check if we need to allocate more memory to 'inf', 'dtime', and 'dsize'
            if (alloc[g.s[j]] == n[g.s[j]].ni) {
                // add 1000 to alloc
                alloc[g.s[j]] += 1000;
                // reallocate memory of 'inf', 'dtime', and 'dsize'
                n[g.s[j]].inf = realloc(n[g.s[j]].inf, alloc[g.s[j]] * sizeof(unsigned int));
                n[g.s[j]].dtime = realloc(n[g.s[j]].dtime, alloc[g.s[j]] * sizeof(unsigned int));
                n[g.s[j]].dsize = realloc(n[g.s[j]].dsize, alloc[g.s[j]] * sizeof(unsigned int));
            }
            
            // compute penalty reduction for size of outbreak
            n[g.s[j]].dsize[n[g.s[j]].ni - 1] = g.ns - n[g.s[j]].dsize[n[g.s[j]].ni - 1];
            
            // set heap and time back to NONE
            n[g.s[j]].heap = n[g.s[j]].time = NONE;
        
        }
        
        // store outbreak size
        g.outbreak_sizes[i] = g.ns;
        
        // increment i
        i++;
	
    }
    
    // since not all simulation runs infect every node,
    // we reallocate memory correctly
    for (i = 0; i < g.n; i++) {
        n[i].inf = realloc(n[i].inf, n[i].ni * sizeof(unsigned int));
        n[i].dtime = realloc(n[i].dtime, n[i].ni * sizeof(unsigned int));
        n[i].dsize = realloc(n[i].dsize, n[i].ni * sizeof(unsigned int));
    }
    
    // sum up all outbreak sizes
    for (i = 0; i < neval; i++) s += (double) g.outbreak_sizes[i];
    
    // print average outbreak size to command line
    printf("Average outbreak size in testing phase is %f\n", s /= neval);
    
    // compute min and max outbreak size
    compute_min(g.outbreak_sizes, neval);
    compute_max(g.outbreak_sizes, neval);
    
    // free memory allocated to g.s and alloc
	free(g.s); free(alloc);

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -