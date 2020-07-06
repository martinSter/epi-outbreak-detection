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
// function to print progress bar

void progress_bar (char label[], int step, int total) {
    
    //progress width
    const int pwidth = 72;

    //minus label len
    int width = pwidth - strlen(label);
    int pos = (step * width) / total;
    
    // compute percentage
    int percent = ( step * 100 ) / total;

    // print label of progress bar
    printf("%s[", label);

    //fill progress bar with =
    for ( int i = 0; i < pos; i++ )  printf( "%c", '=' );

    //fill progress bar with spaces
    printf("%*c", width - pos + 1, ']');
    printf(" %3d%%\r", percent);
    
    // force print operation to complete before the program moves to next task
    fflush(stdout);

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// sort nodes by degree

void sort_by_degree () {
    
    // declare integers
    unsigned int i, j, k, temp1, temp2;
    
    // initialize array on the stack that holds degrees of nodes
    unsigned int degrees[g.n];
    
    // allocate memory to g.deg
    g.deg = malloc(g.n * sizeof(unsigned int));
    
    // loop over all nodes and fill arrays
    for (i = 0; i < g.n; i++) {
        
        // add node i to g.deg
        g.deg[i] = i;
        
        // initialize degree of node i to 0
        degrees[i] = 0;
        
        // loop over all neighbors of node i
        for (j = 0; j < n[i].deg; j++) {
            
            // loop over all contacts with a given neighbor
            for (k = 0; k < n[i].nc[j]; k++) {
                
                // test if some contact with given neighbor is within time period
                if ((n[i].t[j][k] >= g.t_start) && (n[i].t[j][k] <= g.t_end)) {
                    
                    // increment degree count
                    degrees[i]++;
                    
                    // go to next neighbor
                    break;
                    
                }
                
            }
        
        }

    }
 
    // loop over array elements
    for (i = 0; i < g.n; i++) {
        
        // loop over all array elements after i-th element
        for (j = i + 1; j < g.n; j++) {
            
            // swap elements if i-th element is smaller than j-th element
            if (degrees[i] < degrees[j]) {
                
                // store node and degree at i in temp
                temp1 = g.deg[i];
                temp2 = degrees[i];
                
                // set i to j
                g.deg[i] = g.deg[j];
                degrees[i] = degrees[j];
                
                // set j from temp
                g.deg[j] = temp1;
                degrees[j] = temp2;
            }
        }
    }
    
    // print three top nodes
    printf("Node %u has degree %u\n", g.deg[0], degrees[0]);
    printf("Node %u has degree %u\n", g.deg[1], degrees[1]);
    printf("Node %u has degree %u\n", g.deg[2], degrees[2]);
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// sort nodes by number of links

void sort_by_links () {
    
    // declare integers
    unsigned int i, j, k, temp1, temp2;
    
    // initialize array on the stack that holds degrees of nodes
    unsigned int links[g.n];
    
    // allocate memory to g.lin
    g.lin = malloc(g.n * sizeof(unsigned int));
    
    // loop over all nodes and fill arrays
    for (i = 0; i < g.n; i++) {
        
        // add node i to g.lin
        g.lin[i] = i;
        
        // initialize degree of node i to 0
        links[i] = 0;
        
        // loop over all neighbors of node i
        for (j = 0; j < n[i].deg; j++) {
            
            // loop over all contacts with a given neighbor
            for (k = 0; k < n[i].nc[j]; k++) {
                
                // test if some contact with given neighbor is within time period
                if ((n[i].t[j][k] >= g.t_start) && (n[i].t[j][k] <= g.t_end)) {
                    
                    // increment degree count
                    links[i]++;
                }
            }
        }
    }
 
    // loop over array elements
    for (i = 0; i < g.n; i++) {
        
        // loop over all array elements after i-th element
        for (j = i + 1; j < g.n; j++) {
            
            // swap elements if i-th element is smaller than j-th element
            if (links[i] < links[j]) {
                
                // store node and degree at i in temp
                temp1 = g.lin[i];
                temp2 = links[i];
                
                // set i to j
                g.lin[i] = g.lin[j];
                links[i] = links[j];
                
                // set j from temp
                g.lin[j] = temp1;
                links[j] = temp2;
            }
        }
    }
    
    // print three top nodes
    printf("Node %u has %u links\n", g.lin[0], links[0]);
    printf("Node %u has %u links\n", g.lin[1], links[1]);
    printf("Node %u has %u links\n", g.lin[2], links[2]);
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// randomize node order

void shuffle_nodes () {
    
    // declare integers
    unsigned int i, j, t;
    
    // allocate memory to g.ran
    g.ran = malloc(g.n * sizeof(unsigned int));
    
    // initialize g.ran
    for (i = 0; i < g.n; i++) g.ran[i] = i;
    
    // shuffle
    for (i = 0; i < g.n - 1; i++) {
        
        // get random index after i
        j = i + pcg_32_bounded(g.n - i);
        
        // store node at j in t
        t = g.ran[j];
        
        // set node at j to node at i
        g.ran[j] = g.ran[i];
        
        // set node at i to node previously at j
        g.ran[i] = t;
        
    }
   
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// quicksort routines

// utility function to swap two elements 
void swap(int i, int j, unsigned int *out_sizes) {
    
    // declare unsigned integers
    unsigned int temp = out_sizes[i];

    // get values of j at i
    out_sizes[i] = out_sizes[j];
    
    // get values of i at j
    out_sizes[j] = temp;

}
  
/* This function takes last element as pivot, places 
   the pivot element at its correct position in sorted 
    array, and places all smaller (smaller than pivot) 
   to left of pivot and all greater elements to right 
   of pivot */
int partition (int low, int high, unsigned int *out_sizes) {
    
    // set pivot
    unsigned int pivot = out_sizes[high];
    
    // index of smaller element
    int i = (low - 1); 

    // loop from low to high
    for (int j = low; j <= high-1; j++) {
        
        // if current element is smaller than the pivot 
        if (out_sizes[j] < pivot) {
            
            // increment index of smaller element
            i++;
            
            // swap elements i and j
            swap(i, j, out_sizes);
            
        }
        
    }
    
    // swap elements i+1 and high (put pivot in right position)
    swap(i+1, high, out_sizes);
    
    // return partitioning index
    return (i + 1);
    
} 
  
// main function for quickSort
void quickSort(int low, int high, unsigned int *out_sizes) {
    
    // only continue if low is smaller than high
    if (low < high) {
        
        // pi is partitioning index, mg[p] is now at right place
        int pi = partition(low, high, out_sizes);
  
        // separately sort elements before partition and after partition 
        quickSort(low, pi - 1, out_sizes); 
        quickSort(pi + 1, high, out_sizes); 
    
    }
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// median

void compute_median (unsigned int *out_sizes, unsigned int size_array) {
    
    // declare double
    double med;
    
    // quicksort the array
    quickSort(0, size_array - 1, out_sizes); 
  
    // check for even case 
    if (size_array % 2 != 0) med = (double) out_sizes[size_array/2]; 
    else med = (double)(out_sizes[(size_array-1)/2] + out_sizes[size_array/2]) / 2.0;
    
    // print to command line
    printf("Median outbreak size is %f\n", med);
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// max and min

void compute_max (unsigned int *out_sizes, unsigned int size_array) {
    
    // Initialize maximum element 
    unsigned int i, max_val = out_sizes[0];
  
    // Traverse array elements from second and compare every element with current max   
    for (i = 1; i < size_array; i++) if (out_sizes[i] > max_val) max_val = out_sizes[i];
    
    // print to command line
    printf("Max. outbreak size is %u\n", max_val);
    
}

void compute_min (unsigned int *out_sizes, unsigned int size_array) {
    
    // Initialize maximum element 
    unsigned int i, min_val = out_sizes[0];
  
    // Traverse array elements from second and compare every element with current max   
    for (i = 1; i < size_array; i++) if (out_sizes[i] < min_val) min_val = out_sizes[i];
    
    // print to command line
    printf("Min. outbreak size is %u\n", min_val);
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -