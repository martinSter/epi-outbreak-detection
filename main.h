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
#define NSIM 50000

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
	// OUTBREAK STATS
	unsigned int ns, *s;
    // array to store optimal node order
    unsigned int *on;
} GLOBALS;

// struct to define nodes
typedef struct NODE {
	unsigned int deg, *nb; // degree, neighbors
	unsigned int *nc, **t; // ordered number of / list of contact times for bisection search
	unsigned int heap, time; // time is 1st the time of infection (for sorting the heap), then the time of recovery (to check if the node is I or R)
    unsigned int ni, *inf; // ni corresponds to total number of simulations that infect node, inf is an array that contains the indices of simulations that infect node
} NODE;

// struct to define node-marginal-gain pairs
typedef struct MARGINALGAIN {
    unsigned int node;
    unsigned int gain;
} MARGINALGAIN;

// FUNCTION PROTOTYPES FROM ALL SOURCE FILES
// the 'extern' statement is not absolutely necessary (the compiler implicitly assumes it for functions)

// tsir.c
extern void simulate ();

// greedy.c
extern void greedy_max ();

// heap.c
extern void up_heap (unsigned int);
extern void del_root ();

// misc.c
extern void init_rng ();
extern void read_data (FILE *);
extern unsigned int exptime ();
extern void progress_bar (char label[], int step, int total);

// quick.c
extern void quick (unsigned int);

// pcg_rnd.c
extern uint16_t pcg_16 ();
extern uint32_t pcg_32 ();
extern uint32_t pcg_32_bounded ();
extern void pcg_init ();

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
