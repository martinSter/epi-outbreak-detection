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
// ...

void greedy_max () {
    
    // declare unsigned int i
    unsigned int i, idx = 0, last = g.n - 1;
    
    // allocate memory to g.on (array storing the optimal node order)
    g.on = malloc(g.n * sizeof(unsigned int));
    
    // initialize number of detected scenarios to 0
    g.nd = 0;
    
    // allocate memory to g.detected
    g.detected = malloc(NSIM * sizeof(unsigned int));
    
    // allocate memory to mg
    mg = calloc(g.n, sizeof(MARGINALGAIN));
    
    // initialize mg
    for (i = 0; i < g.n; i++) {
        mg[i].node = i;
        mg[i].gain = n[i].ni;
    }
    
    // sort mg (ascending order)
    // first argument is first index and second argument is last index
    quickSort(0, last);
    
    printf("Node %u has max. marginal gain %u\n", mg[g.n-1].node, mg[g.n-1].gain);
    printf("Node %u has max. marginal gain %u\n", mg[g.n-2].node, mg[g.n-2].gain);
    printf("Node %u has max. marginal gain %u\n", mg[g.n-3].node, mg[g.n-3].gain);
    
    while (last >= 0) {
        
        // special procedure for first node
        if (idx == 0) {
            // add node with largest marginal gain to g.on
            g.on[idx++] = mg[last].node;
            // save scenarios detected by this node
            // ...
            // decrease last by 1
            last--;
        }
        
        break;
        
        // get largest node
        
        // recompute its marginal gain
        
        
        // resort mg
        
        // if node remains at the top of array
        
            // add it to g.on
            
            // decrease last by 1
            
            // save scenarios detected by this node
            
            // empty active
        
        // else
            
            // if node is not in active, add it
            
            // else, randomly sample node from all active nodes
            
            // etc.
    
    
    
    
    }
    
    
    
    
    // deallocate memory
    free(mg);
    
}









// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
