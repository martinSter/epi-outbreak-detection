// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for greedy maximization

// include header file
#include "main.h"

// declare g as a GLOBALS struct
extern GLOBALS g;

// declare n as a array of NODE structs
extern NODE *n;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// binary-heap routines (max-heap)

// up-heap operation adjusted from heap.c (max-heap instead of min-heap)
void heap_up (unsigned int start) {
    
    // initialized integers
	unsigned int above, here = start, mem = g.heap[start];

    // while loop until here is at index 1
	while (here > 1) {
        
        // find parent of here
		above = here >> 1; // = here / 2 (to find parent)
		
        // break loop if marginal gain of node is smaller than the one of parent
		if (n[mem].mg <= n[g.heap[above]].mg) break;
        
        // swap child and parent
		g.heap[here] = g.heap[above];
        
        // set heap of parent to here
		n[g.heap[here]].heap = here;
		
        // now, here is above and continue loop
		here = above;
	}
	
    // when loop is done, modify heap of node
	n[g.heap[here] = mem].heap = here;
    
}

// down-heap operation adjusted from heap.c
void heap_down (unsigned int here) {
    
    // initialize unsigned integers
	unsigned int utmp, largest = here;
    
    // find left and right child
	unsigned int left = here << 1; // = here * 2
	unsigned int right = left | 1; // = left + 1

	// if the heap property is violated vs the children, find the largest child 
	if ((left <= g.nheap) && (n[g.heap[left]].mg > n[g.heap[largest]].mg)) largest = left;
	if ((right <= g.nheap) && (n[g.heap[right]].mg > n[g.heap[largest]].mg)) largest = right;

    // if largest is still here, the heap property is not violated and we leave function
	if (largest == here) return;

	// swap largest and here in heap
	utmp = g.heap[largest];
	g.heap[largest] = g.heap[here];
	g.heap[here] = utmp;

    // modify heap positions of nodes
	n[g.heap[largest]].heap = largest;
	n[g.heap[here]].heap = here;

	// continue checking below
	heap_down(largest);
}

// extract operation adjusted from heap.c
void remove_root () {

    // set heap position of root to END
	n[g.heap[1]].heap = END;
    
    // set last element in heap as new root and decrement g.nheap
	g.heap[1] = g.heap[g.nheap--];
    
    // set heap position of new root to 1
	n[g.heap[1]].heap = 1;
    
    // down-heap new root
	heap_down(1);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// utility functions

// store the detected scenarios in g.detected
void store_detected_dl (unsigned int v) {
    
    // declare unsigned int i
    unsigned int i;
    
    // set all scenarios that node v detects to 1
    for (i = 0; i < n[v].ni; i++) g.detected[n[v].inf[i]] = 1;
    
}

// recompute the marginal gain of node v
void recompute_mg_dl (unsigned int v) {
    
    // declare unsigned int i
    unsigned int i;
    
    // set current marginal gain of node v to 0
    n[v].mg = 0;
    
    // recompute marginal gain of node v
    for (i = 0; i < n[v].ni; i++) if (g.detected[n[v].inf[i]] == 0) n[v].mg++;
    
}

// store the detected scenarios in g.detected
void store_detected_dt (unsigned int v) {
    
    // declare unsigned int i
    unsigned int i;
    
    // set all scenarios that node v detects to 1
    for (i = 0; i < n[v].ni; i++) {
        
        // if penalty reduction (pr) of v is larger than current pr, then replace it
        if (n[v].dtime[i] > g.detected[n[v].inf[i]]) g.detected[n[v].inf[i]] = n[v].dtime[i];        
        
    }
    
}

// recompute the marginal gain of node v
void recompute_mg_dt (unsigned int v) {
    
    // declare unsigned int i
    unsigned int i;
    
    // set current marginal gain of node v to 0
    n[v].mg = 0;
    
    // recompute marginal gain of node v
    for (i = 0; i < n[v].ni; i++) {
        
        // compute difference between pa of node v and previous pa for scenario i
        int temp = n[v].dtime[i] - g.detected[n[v].inf[i]];
        
        // check if penalty reduction of node v is larger than current penalty reduction (for some scenario)
        if (temp > 0) {
            
            // increment marginal gain of node v by difference
            n[v].mg += temp;
            
        }
        
    }
        
}

// store the detected scenarios in g.detected
void store_detected_pa (unsigned int v) {
    
    // declare unsigned int i
    unsigned int i;
    
    // set all scenarios that node v detects to 1
    for (i = 0; i < n[v].ni; i++) {
        
        // if penalty reduction (pr) of v is larger than current pr, then replace it
        if (n[v].dsize[i] > g.detected[n[v].inf[i]]) g.detected[n[v].inf[i]] = n[v].dsize[i];        
        
    }
    
}

// recompute the marginal gain of node v
void recompute_mg_pa (unsigned int v) {
    
    // declare unsigned int i
    unsigned int i;
    
    // set current marginal gain of node v to 0
    n[v].mg = 0;
    
    // recompute marginal gain of node v
    for (i = 0; i < n[v].ni; i++) {
        
        // compute difference between pa of node v and previous pa for scenario i
        int temp = n[v].dsize[i] - g.detected[n[v].inf[i]];
        
        // check if penalty reduction of node v is larger than current penalty reduction (for some scenario)
        if (temp > 0) {
            
            // increment marginal gain of node v by difference
            n[v].mg += temp;
            
        }
        
    }
        
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// greedy max. (main functions)

// ************************
// detection likelihood (dl)
void greedy_max_dl () {
    
    // declare unsigned int i
    unsigned int i, idx = 0, cn, sum = 0;
    
    // allocate memory to g.on (array storing the optimal node order)
    g.on = malloc(g.n * sizeof(unsigned int));
    
    // allocate memory to g.detected
    g.detected = calloc(NSIM, sizeof(unsigned int));
    
    // initialize marginal gains to number of infected scenarios for every node
    for (i = 0; i < g.n; i++) n[i].mg = n[i].ni;
    
    // set nheap to 0
    g.nheap = 0;
    
    // initialize max-heap
    for (i = 0; i < g.n; i++) {
        
        // add node i to heap
        g.heap[++g.nheap] = i;
        
        // modify heap position of node i
        n[i].heap = g.nheap;
        
        // up-heap node i
        heap_up(n[i].heap);
        
    }
    
    // add first node in heap to g.on
    g.on[idx++] = g.heap[1];
    
    // print marginal gain of node added to g.on
    printf("Node %u has been added with marginal gain = %lu\n", g.heap[1], n[g.heap[1]].mg);
    
    // store scenarios detected by root
    store_detected_dl(g.heap[1]);
    
    // remove root from heap
    remove_root();
    
    // loop as long as heap contains nodes
    while (g.nheap > 0) {
        
        // store the current root node
        cn = g.heap[1];
        
        // recompute marginal gain of root node
        recompute_mg_dl(g.heap[1]);
        
        // down-heap root if necessary
        heap_down(1);
        
        // if former root is still on top of heap, we add it to g.on
        if (g.heap[1] == cn) {
            
            // add root to g.on
            g.on[idx++] = g.heap[1];
            
            // store scenarios detected by root
            store_detected_dl(g.heap[1]);
            
            // remove root from heap
            remove_root();
            
        }
        
    }
    
    // make sure we detect all scenarios
    for (i = 0; i < NSIM; i++) sum += g.detected[i];
    printf("Number of detected scenarios is %u\n", sum);
    
    // free memory allocated to g.detected
    free(g.detected);
    
}

// ************************
// detection time (dt)
void greedy_max_dt () {
    
    // declare unsigned int i
    unsigned int i, j, idx = 0, cn;
    
    // allocate memory to g.on (array storing the optimal node order)
    g.dt = malloc(g.n * sizeof(unsigned int));
    
    // allocate memory to g.detected
    g.detected = calloc(NSIM, sizeof(unsigned int));
    
    // initialize marginal gains of nodes
    for (i = 0; i < g.n; i++) {
        
        // set marginal gains to 0
        n[i].mg = 0;
        
        // add up all penalty reductions
        for (j = 0; j < n[i].ni; j++) n[i].mg += n[i].dtime[j];
        
    }
    
    // set nheap to 0
    g.nheap = 0;
    
    // initialize max-heap
    for (i = 0; i < g.n; i++) {
        
        // add node i to heap
        g.heap[++g.nheap] = i;
        
        // modify heap position of node i
        n[i].heap = g.nheap;
        
        // up-heap node i
        heap_up(n[i].heap);
        
    }
    
    // add first node in heap to g.dt
    g.dt[idx++] = g.heap[1];
    
    // print marginal gain of node added to g.on
    printf("Node %u has been added with marginal gain = %lu\n", g.heap[1], n[g.heap[1]].mg);
    
    // store scenarios detected by root
    store_detected_dt(g.heap[1]);
    
    // remove root from heap
    remove_root();
    
    // loop as long as heap contains nodes
    while (g.nheap > 0) {
        
        // store the current root node
        cn = g.heap[1];
        
        // recompute marginal gain of root node
        recompute_mg_dt(g.heap[1]);
        
        // down-heap root if necessary
        heap_down(1);
        
        // if former root is still on top of heap, we add it to g.on
        if (g.heap[1] == cn) {
            
            // add root to g.on
            g.dt[idx++] = g.heap[1];
            
            // store scenarios detected by root
            store_detected_dt(g.heap[1]);
            
            // remove root from heap
            remove_root();
            
        }
        
    }
    
    // free memory allocated to g.detected
    free(g.detected);
    
}

// ************************
// population affected (pa)
void greedy_max_pa () {
    
    // declare unsigned int i
    unsigned int i, j, idx = 0, cn;
    
    // allocate memory to g.on (array storing the optimal node order)
    g.pa = malloc(g.n * sizeof(unsigned int));
    
    // allocate memory to g.detected
    g.detected = calloc(NSIM, sizeof(unsigned int));
    
    // initialize marginal gains of nodes
    for (i = 0; i < g.n; i++) {
        
        // set marginal gains to 0
        n[i].mg = 0;
        
        // add up all penalty reductions
        for (j = 0; j < n[i].ni; j++) n[i].mg += n[i].dsize[j];
        
    }
    
    // set nheap to 0
    g.nheap = 0;
    
    // initialize max-heap
    for (i = 0; i < g.n; i++) {
        
        // add node i to heap
        g.heap[++g.nheap] = i;
        
        // modify heap position of node i
        n[i].heap = g.nheap;
        
        // up-heap node i
        heap_up(n[i].heap);
        
    }
    
    // add first node in heap to g.dt
    g.pa[idx++] = g.heap[1];
    
    // print marginal gain of node added to g.on
    printf("Node %u has been added with marginal gain = %lu\n", g.heap[1], n[g.heap[1]].mg);
    
    // store scenarios detected by root
    store_detected_pa(g.heap[1]);
    
    // remove root from heap
    remove_root();
    
    // loop as long as heap contains nodes
    while (g.nheap > 0) {
        
        // store the current root node
        cn = g.heap[1];
        
        // recompute marginal gain of root node
        recompute_mg_pa(g.heap[1]);
        
        // down-heap root if necessary
        heap_down(1);
        
        // if former root is still on top of heap, we add it to g.on
        if (g.heap[1] == cn) {
            
            // add root to g.on
            g.pa[idx++] = g.heap[1];
            
            // store scenarios detected by root
            store_detected_pa(g.heap[1]);
            
            // remove root from heap
            remove_root();
            
        }
        
    }
    
    // free memory allocated to g.detected
    free(g.detected);
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
