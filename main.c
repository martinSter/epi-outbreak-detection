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
	unsigned int i, j, k;
    
    // declare d, initialize s1 and s2 (for averages)
	double d;
    
    // declare pointer to FILE
	FILE *fp;
	
	// just a help message
	if ((argc < 5) || (argc > 5)) {
		fprintf(stderr, "usage: ./main [data file] [beta] [nu (units of the duration of the data)] [seed]\n");
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

	// compute 1/log(1-beta) --> inverse probability transform
	d = 1.0 / log(1.0 - atof(argv[2]));
	
    // create a large number (65537) of random numbers created according to geometric distribution
    for (i = 0; i < 0x10000; i++)
		g.rnd2inx[i] = (unsigned short) floor(d * log((i + 1) / 65536.0));
	
    // compute recovery parameter
    g.recovery_scale = g.dur / atof(argv[3]);
    
    // initialize random number generator
    g.state = (uint64_t) strtoull(argv[4], NULL, 10);

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
    
    for (i = 0; i < n[15].ni; i++) printf("Simulation run: %u, time penalty reduction: %u, size penalty reduction: %u\n", n[15].inf[i], n[15].dtime[i], n[15].dsize[i]);
    
    
    // for (i = 1; i < 10; i++) printf("Sim. ID: %d\n", n[11134].inf[i]);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // GREEDY MAXIMIZATION
    
    greedy_max_dl();
    greedy_max_dt();
    greedy_max_pa();
    
    // free memory allocated to inf, dtime, and dsize
    for (i = 0; i < g.n; i++) {
        free(n[i].inf);
        free(n[i].dtime);
        free(n[i].dsize);
    }
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // SORT NODES BY DEGREE
    
    sort_by_degree();
    
    printf("Node %u has degree %u\n", g.deg[0], n[g.deg[0]].deg);
    printf("Node %u has degree %u\n", g.deg[1], n[g.deg[1]].deg);
    printf("Node %u has degree %u\n", g.deg[2], n[g.deg[2]].deg);
    
    printf("Node %u has degree %u\n", 11134, n[11134].deg);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // RANDOMLY SHUFFLE NODES
    
    shuffle_nodes();
    
    printf("Random node %u\n", g.ran[0]);
    printf("Random node %u\n", g.ran[1]);
    printf("Random node %u\n", g.ran[2]);
    
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // EVALUATION
        
    // size of evaluation set
    unsigned int neval = 10000;
    
    // allocate memory to g.res_greedy
    g.res_greedy = calloc(g.n, sizeof(unsigned int));
    g.res_degree = calloc(g.n, sizeof(unsigned int));
    g.res_random = calloc(g.n, sizeof(unsigned int));
    
    // create evaluation set of outbreak scenarios
    simulate_eval(neval);
    
    // allocate memory to g.detected
    g.detected = calloc(neval, sizeof(unsigned int));
    
    // RESULTS GREEDY
    
    // loop over g.on and find marginal improvements
    for (i = 0; i < g.n; i++) {
        
        // set all scenarios that node v detects to 1
        for (j = 0; j < n[g.on[i]].ni; j++) g.detected[n[g.on[i]].inf[j]] = 1;
        
        // sum number of detected cases
        for (k = 0; k < neval; k++) g.res_greedy[i] += g.detected[k]; 
        
    }
    
    // RESULTS DEGREE
    
    // set all elements in g.detected back to 0
    memset(g.detected, 0, neval*sizeof(unsigned int));
    
    // loop over g.on and find marginal improvements
    for (i = 0; i < g.n; i++) {
        
        // set all scenarios that node v detects to 1
        for (j = 0; j < n[g.deg[i]].ni; j++) g.detected[n[g.deg[i]].inf[j]] = 1;
        
        // sum number of detected cases
        for (k = 0; k < neval; k++) g.res_degree[i] += g.detected[k]; 
        
    }
    
    // RESULTS RANDOM
    
    // set all elements in g.detected back to 0
    memset(g.detected, 0, neval*sizeof(unsigned int));
    
    // loop over g.on and find marginal improvements
    for (i = 0; i < g.n; i++) {
        
        // set all scenarios that node v detects to 1
        for (j = 0; j < n[g.ran[i]].ni; j++) g.detected[n[g.ran[i]].inf[j]] = 1;
        
        // sum number of detected cases
        for (k = 0; k < neval; k++) g.res_random[i] += g.detected[k]; 
        
    }
    
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // EXPORT RESULTS
    
    // open network data file
	fp = fopen("nodes.txt", "w");
	if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return 1;
	}
    
    // print optimal nodes to file
    for (i = 0; i < g.n; i++) fprintf(fp, "%u;%u;%u;%u;%u\n", g.on[i], g.dt[i], g.pa[i], g.deg[i], g.ran[i]);
    
    // print data
    // for (i = 0; i < g.n; i++) fprintf(fp, "%u;%u;%u;%u;%u;%u\n", g.on[i], g.res_greedy[i], g.deg[i], g.res_degree[i], g.ran[i], g.res_random[i]);
    
    // close file
	fclose(fp);
    
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
        free(n[i].dtime);
        free(n[i].dsize);
	}
    
    // free array n of NODE structs and heap and s (only heap and s are defined as pointers in GLOBALS)
	free(n); free(g.heap); free(g.detected);
    free(g.on); free(g.dt); free(g.pa); free(g.deg); free(g.ran);
    free(g.res_greedy); free(g.res_degree); free(g.res_random);
	 
	return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
