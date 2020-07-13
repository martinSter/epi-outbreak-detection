// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for temporal network SIR by Petter Holme (2018)

// include external libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>

// number of simulation runs
#define NSIM 5000000

// minimal outbreak size for scenario generation
#define MIN_OUTSIZE 2

// NONE and END are used in various ways, the only purpose of NONE < END is for the S(x) macro
// UINT_MAX is the maximum value for an object of type unsigned int.
// Those integers are effectively used to assign time infinity.
#define NONE (UINT_MAX - 1)
#define END UINT_MAX

// is x susceptible?
#define S(x) (n[(x)].heap < END)

// auxiliary macro to compute square of x
#define SQ(x) ((x) * (x))

// struct to define global parameters
typedef struct GLOBALS {
	// INPUT PARAMETERS
	double recovery_scale; // recovery time scale, auxiliary value for infection probs
	unsigned short rnd2inx[0x10001]; // mapping 16-bit random number to index
	// NETWORK SPECS (number of nodes, duration)
	unsigned int n, dur;
	// OTHER GLOBALS
	unsigned int nheap, *heap;
	// FOR RND (random number generation)
	uint64_t state;
	uint32_t rmem;
	unsigned int cutoff_source, cutoff_dur; // to get the probabilities right . .
    // START AND END TIME
    unsigned int t_start, t_end;
	// OUTBREAK STATS
	unsigned int ns, sim_id, *s;
    // STORE ALL OUTBREAK SIZES
    unsigned int *outbreak_sizes;
    // array to store optimal node order
    unsigned int *on, *dt, *pa, *deg, *lin, *ran;
    // store number of detected scenarios and scenario IDs
    unsigned int nd, *detected;
    // store results
    unsigned int *res_greedy, *res_greedy_dt, *res_greedy_pa;
    unsigned int *res_degree, *res_degree_dt, *res_degree_pa;
    unsigned int *res_links, *res_links_dt, *res_links_pa;
    unsigned int *res_random, *res_random_dt, *res_random_pa;
} GLOBALS;

// struct to define nodes
typedef struct NODE {
	unsigned int deg, *nb; // degree, neighbors
	unsigned int *nc, **t; // ordered number of / list of contact times for bisection search
	unsigned int heap, time; // time is 1st the time of infection (for sorting the heap), then the time of recovery (to check if the node is I or R)
    unsigned int ni, *inf; // ni corresponds to total number of simulations that infect node, mg is the current marginal gain, inf is an array that contains the indices of simulations that infect node
    unsigned int *dtime, *dsize; // arrays to store the detection time and the size of the outbreak at detection time
    unsigned long int mg; // marginal gain of a node
} NODE;

// struct to define node-marginal-gain pairs
typedef struct MARGINALGAIN {
    unsigned int node;
    unsigned int gain;
} MARGINALGAIN;

// FUNCTION PROTOTYPES FROM ALL SOURCE FILES
// the 'extern' statement is not absolutely necessary (the compiler implicitly assumes it for functions)

// tsir.c
extern void simulate (unsigned int t_start, unsigned int t_end);
extern void simulate_eval (unsigned int neval, unsigned int t_start, unsigned int t_end);

// greedy.c
extern void greedy_max_dl ();
extern void greedy_max_dt ();
extern void greedy_max_pa ();

// heap.c
extern void up_heap (unsigned int);
extern void del_root ();

// misc.c
extern void read_data (FILE *);
extern unsigned int exptime ();
extern void progress_bar (char label[], int step, int total);
extern void sort_by_degree ();
extern void sort_by_links ();
extern void shuffle_nodes ();
extern void quickSort(int low, int high, unsigned int *out_sizes);
extern void compute_median (unsigned int *out_sizes, unsigned int size_array);
extern void compute_min (unsigned int *out_sizes, unsigned int size_array);
extern void compute_max (unsigned int *out_sizes, unsigned int size_array);

// quick.c
extern void quick (unsigned int);

// pcg_rnd.c
extern uint16_t pcg_16 ();
extern uint32_t pcg_32 ();
extern uint32_t pcg_32_bounded (uint32_t bound);
extern uint32_t pcg_32_bounded_ul (uint32_t lower, uint32_t upper);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
