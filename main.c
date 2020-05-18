// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for SIR on temporal networks by Petter Holme (2018)

// include header file
#include "main.h"

// declare g as a GLOBALS struct
GLOBALS g;

// declare n as a array of NODE structs
NODE *n;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main function handling i/o

int main (int argc, char *argv[]) {
    
    // declare integers i and j
	unsigned int i, j;
    
    // declare d, initialize s1 and s2 (for averages)
	double d;
    
    // declare pointer to FILE
	FILE *fp;
	
	// just a help message
	if ((argc < 4) || (argc > 4)) {
		fprintf(stderr, "usage: ./main [data file] [beta] [nu (units of the duration of the data)]\n");
		return 1;
	}
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// open network data file
	fp = fopen(argv[1], "r");
	if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open '%s'\n", argv[1]);
		return 1;
	}
    
    // call read_data function
	read_data(fp);
    
    // close file
	fclose(fp);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // initialize random number generator
	pcg_init();

	// compute 1/log(1-beta) --> inverse probability transform
	d = 1.0 / log(1.0 - atof(argv[2]));
	
    // create a large number (65537) of random numbers created according to geometric distribution
    for (i = 0; i < 0x10000; i++)
		g.rnd2inx[i] = (unsigned short) floor(d * log((i + 1) / 65536.0));
	
    // compute recovery parameter
    g.recovery_scale = g.dur / atof(argv[3]);

	// allocating the heap (N + 1) because its indices are 1,...,N
    // g.heap is a pointer, so here we assign it the address where it will point to
	g.heap = malloc((g.n + 1) * sizeof(unsigned int));
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // SIMULATION OF OUTBREAKS
    
    // simulate NSIM times and store scenario ID's in 'n'
    simulate();
    
    // Initialize maximum element 
    unsigned int max_val = n[0].ni, max_node = 0; 
  
    // Traverse array elements from second and 
    // compare every element with current max   
    for (i = 1; i < g.n; i++) {
        if (n[i].ni > max_val) {
            max_val = n[i].ni;
            max_node = i;
        }
    }
    
    printf("\nNode %u discovers the max. number of scenarios (%u)\n", max_node, max_val);
    
    for (i = 1; i < 10; i++) printf("Sim. ID: %d\n", n[11134].inf[i]);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // GREEDY MAXIMIZATION
    
    greedy_max();
    
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // CLEANING UP
    
	// deallocate memory
	for (i = 0; i < g.n; i++) {
		for (j = 0; j < n[i].deg; j++) free(n[i].t[j]);
        // free all arrays in NODE struct n
		free(n[i].nb);
		free(n[i].nc);
		free(n[i].t);
        free(n[i].inf);
	}
    
    // free array n of NODE structs and heap and s (only heap and s are defined as pointers in GLOBALS)
	free(n); free(g.heap); free(g.on);
	 
	return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
