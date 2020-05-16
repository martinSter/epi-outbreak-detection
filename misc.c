// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for temporal network SIR by Petter Holme (2018)
// miscellaneous routines for tsir

// include header file
#include "main.h"

// declare external variables (they have been defined in 'tsir.c')
extern GLOBALS g;
extern NODE *n; // n must be an array of NODE structs

// declare and define pointer 'alloc'
unsigned int *alloc; // this must also be an array of integers

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// giving exponential random numbers with mean 'scale'

unsigned int exptime () {
    // this creates a uniform random number d between 0 and 1
	double d = (4294967296 - pcg_32_bounded(4294967295)) / 4294967296.0;
    // use inverse probability transform of exponential to randomly sample and take floor, cast as unsigned int
	return (unsigned int) floor(-(log(d) * g.recovery_scale));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// gets the index of you in me's adjacency list

unsigned int get_index (unsigned int me, unsigned int you) {
	
    // define int i
    unsigned int i;

    // loop over neighbors of 'me' and find index i where neighbor is 'you'
    // if we find neighbor, we return its index
	for (i = 0; i < n[me].deg; i++) if (n[me].nb[i] == you) return i;

	// this is procedure if neighbor is not yet in 'nb'
    // note: alloc is initialized with all zeros
    // note: we allocate more space than we need (double deg) so that we
    // do not have to allocate more space in every step
	if (alloc[me] <= n[me].deg) {
		// get value of alloc at index 'me'
        i = alloc[me];
        // set the value at alloc[me] to double the degree if degree is positive or 
        // to 1 otherwise (i.e. if deg == 0)
		alloc[me] = (n[me].deg > 0) ? 2 * n[me].deg : 1;
        // increase the allocated space to arrays with the value in alloc[me]
		n[me].nb = realloc(n[me].nb, alloc[me] * sizeof(unsigned int));
		n[me].nc = realloc(n[me].nc, alloc[me] * sizeof(unsigned int));
        // since t is a pointer to a pointer, we use the size of a pointer to int
        // effectively, t is an array of arrays containing the contact times with each neighbor
		n[me].t = realloc(n[me].t, alloc[me] * sizeof(unsigned int *));
		// for loop that does not start at 0 but at previous value of alloc[me], see above
        for ( ; i < alloc[me]; i++) {
            // initialize new positions in nc and nb to 0
			n[me].nc[i] = n[me].nb[i] = 0;
            // initialize new arrays in t to NULL
			n[me].t[i] = NULL;
		}
	}

    // add 'you' to 'nb' and increase 'deg' by one
	n[me].nb[n[me].deg++] = you;

    // return index of new neighbor
	return n[me].deg - 1;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// reads the network, assumes a contact list with vertex label 0,N-1,
// ordered in time. If your network has nodes with zero contacts, make sure
// that none of them is the node with largest index

void read_data (FILE *fp) {
    
    // declare integers
	unsigned int i, me, you;

    // initialize number of nodes in g to 0
	g.n = 0;

	// scan the system size
    // loop over lines in fp
    // that is why the node with largest integer should have contacts,
    // if it does not, then the system size will not be collected correctly
	while (2 == fscanf(fp, "%u %u %*u\n", &me, &you)) {
		// set g.n always to largest value
        if (g.n < me) g.n = me;
		if (g.n < you) g.n = you;
	}

    // adds 1 to n since node indices start at 0
	g.n++;

    // initialize array of NODE structs in heap
    // with calloc we initialize integers to 0
	n = calloc(g.n, sizeof(NODE));
    
    // initialize array of size system (number of nodes)
    // with calloc we initialize integers to 0
    // this is used to manage memory allocation of nb, nc, and t
	alloc = calloc(g.n, sizeof(unsigned int));

    // sets the file position to the beginning
	rewind(fp);

	// scan the degrees (here we implicitly assume undirected network)
    // loop over lines in fp
	while (2 == fscanf(fp, "%u %u %*u\n", &me, &you)) {
		// find index of 'you' in adjacency list of 'me'
        // accounting for deg and nb is taken care of in 'get_index()'
        i = get_index(me, you);
		// increase number of contacts with node 'you' by 1
        n[me].nc[i]++;
		// find index of 'me' in adjacency list of 'you'
        // accounting for deg and nb is taken care of in 'get_index()'
        i = get_index(you, me);
        // increase number of contacts with node 'me' by 1
		n[you].nc[i]++;
	}

    // sets the file position to the beginning
	rewind(fp);

    // prepare arrays of contact times
    // loop over int me (over all nodes in network)
	for (me = 0; me < g.n; me++) {
        // loop over all contacts 'i' of node 'me'
		for (i = 0; i < n[me].deg; i++) {
            // allocate memory for array of contact times
			n[me].t[i] = malloc(n[me].nc[i] * sizeof(unsigned int));
            // set number of contacts back to 0
            // this is needed because we use nc to increment in next code block
			n[me].nc[i] = 0;
		}
	}

	// scan the times
    // loop over lines in fp (this time considering contact times)
    // since input edgelist is sorted by date, the last g.dur will be length of period
	while (3 == fscanf(fp, "%u %u %u\n", &me, &you, &g.dur)) {
		// find index of 'you' in adjacency list of 'me'
        i = get_index(me, you);
        // set contact times
		n[me].t[i][n[me].nc[i]++] = g.dur;
        // find index of 'me' in adjacency list of 'you'
		i = get_index(you, me);
        // set contact times
		n[you].t[i][n[you].nc[i]++] = g.dur;
	}

	// allocate adjacency lists
    // this is needed because we might have allocated too much memory in the 'get_index()' part
    // here, we therefore allocate the correct amount of memory
    // loop over all nodes in network
	for (i = 0; i < g.n; i++) {
		n[i].nb = realloc(n[i].nb, n[i].deg * sizeof(unsigned int));
		n[i].nc = realloc(n[i].nc, n[i].deg * sizeof(unsigned int));
		n[i].t = realloc(n[i].t, n[i].deg * sizeof(unsigned int *));
        // for every node i, sort t in decreasing order of its last element
		quick(i);
	}
    
    // initialize ni and inf
    for (i = 0; i < g.n; i++) {
        // set ni (number of times node i gets infected) to 0
        n[i].ni = 0;
        // allocate memory for arrays that store simulation runs that infect node i
        n[i].inf = calloc(NSIM, sizeof(unsigned int));
    }

    // free memory
	free(alloc);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
