// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for greedy maximization

// include header file
#include "main.h"

// declare g as a GLOBALS struct
extern GLOBALS g;

// declare n as a array of NODE structs
extern NODE *n;

// declare array of type MARGINALGAIN
MARGINALGAIN *mg;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// quickSort routines

// utility function to swap two elements 
void swap(int i, int j) {
    
    // declare unsigned integers
    unsigned int temp1 = mg[i].node, temp2 = mg[i].gain;

    // get values of j at i
    mg[i].node = mg[j].node;
    mg[i].gain = mg[j].gain;
    
    // get values of i at j
    mg[j].node = temp1;
    mg[j].gain = temp2;

}
  
/* This function takes last element as pivot, places 
   the pivot element at its correct position in sorted 
    array, and places all smaller (smaller than pivot) 
   to left of pivot and all greater elements to right 
   of pivot */
int partition (int low, int high) {
    
    // set pivot
    unsigned int pivot = mg[high].gain;
    
    // index of smaller element
    int i = (low - 1); 

    // loop from low to high
    for (int j = low; j <= high-1; j++) {
        
        // if current element is smaller than the pivot 
        if (mg[j].gain < pivot) {
            
            // increment index of smaller element
            i++;
            
            // swap elements i and j
            swap(i, j);
            
        }
        
    }
    
    // swap elements i+1 and high (put pivot in right position)
    swap(i+1, high);
    
    // return partitioning index
    return (i + 1);
    
} 
  
// main function for quickSort
void quickSort(int low, int high) {
    
    // only continue if low is smaller than high
    if (low < high) {
        
        // pi is partitioning index, mg[p] is now at right place
        int pi = partition(low, high);
  
        // separately sort elements before partition and after partition 
        quickSort(low, pi - 1); 
        quickSort(pi + 1, high); 
    
    }
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void update_scenarios(unsigned int node) {
    
    // set i to current value of g.nd
    unsigned int i = g.nd;
    
    // update g.nd with new scenarios from 'node'
    g.nd += n[node].ni;
    
    // add 
    for (; i < g.nd; i++) g.detected[i] = n[node].inf[i];
    
    // sort g.detected
    
    
}

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
void store_detected (unsigned int v) {
    
    // declare unsigned int i
    unsigned int i;
    
    // set all scenarios that node v detects to 1
    for (i = 0; i < n[v].ni; i++) g.detected[n[v].inf[i]] = 1;
    
}

// recompute the marginal gain of node v
void recompute_mg (unsigned int v) {
    
    // declare unsigned int i
    unsigned int i;
    
    // set current marginal gain of node v to 0
    n[v].mg = 0;
    
    // recompute marginal gain of node v
    for (i = 0; i < n[v].ni; i++) if (g.detected[n[v].inf[i]] == 0) n[v].mg++;
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// greedy max. (main function)

void greedy_max () {
    
    // declare unsigned int i
    unsigned int i, idx = 0, cn, sum = 0;
    
    // allocate memory to g.on (array storing the optimal node order)
    g.on = malloc(g.n * sizeof(unsigned int));
    
    // initialize number of detected scenarios to 0
    g.nd = 0;
    
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
    printf("Node %u has been added with marginal gain = %u\n", g.heap[1], n[g.heap[1]].mg);
    
    // store scenarios detected by root
    store_detected(g.heap[1]);
    
    // remove root from heap
    remove_root();
    
    // loop as long as heap contains nodes
    while (g.nheap > 0) {
        
        // store the current root node
        cn = g.heap[1];
        
        // recompute marginal gain of root node
        recompute_mg(g.heap[1]);
        
        // down-heap root if necessary
        heap_down(1);
        
        // if former root is still on top of heap, we add it to g.on
        if (g.heap[1] == cn) {
            
            // add root to g.on
            g.on[idx++] = g.heap[1];
            
            // print marginal gain of node added to g.on
            // printf("Node %u has been added with marginal gain = %u\n", g.heap[1], n[g.heap[1]].mg);
            
            // store scenarios detected by root
            store_detected(g.heap[1]);
            
            // remove root from heap
            remove_root();
            
        }
        
    }
    
    for (i = 0; i < NSIM; i++) sum += g.detected[i];
    printf("Number of detected scenarios is %u\n", sum);
    
    printf("Max of heap is %u\n", g.nheap);
    
    // free memory allocated to inf
    for (i = 0; i < g.n; i++) {
        free(n[i].inf);
        free(n[i].dtime);
        free(n[i].dsize);
    }
    
    // free memory allocated to g.detected
    free(g.detected);
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
