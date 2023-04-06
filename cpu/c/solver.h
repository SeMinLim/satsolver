#include <stdio.h>
#include <stdlib.h>

#define ChildLeft(x) (x << 1 | 1)
#define ChildRight(x) ((x + 1) << 1)
#define Parent(x) ((x - 1) >> 1)

#define Value(literal) (literal > 0 ? value[literal] : -value[-literal])
#define WatchedPointers(id) (watchedPointers[vars + id])
#define WatchedPointersSize(id) (watchedPointersSize[vars + id])

#define NumVars 99
#define NumClauses 264
#define MaxNumLits 3


// Function for soting activities
typedef struct GreaterActivity { 
    	const double *activity;     
    	int compare ( int a, int b ) { 
		if ( activity[a] > activity[b] ) return 1; 
		else return 0;
	}
    	void initialize( const double *s ) { 
		activity = s;
	}
} GreaterActivity;


// Heap data structure
typedef struct Heap {
	GreaterActivity g;
	int heap[NumVars+1];
	int heapSize = 0;
	int pos[NumVars+1];
	int posSize = 0;
    
    	void up( int v ) {
        	int x = heap[v];
		int p = Parent(v);
		// Child > Parent -> True
        	while ( v && g.compare(x, heap[p]) ) {
       			heap[v] = heap[p];
			pos[heap[p]] = v;
            		v = p; 
			p = Parent(p);
        	}
        	heap[v] = x;
		pos[x] = v;
    	}

    	void down( int v ) {
        	int x = heap[v];
        	while ( ChildLeft(v) < heapSize ){
            		// Pick the bigger one among left and right child
			int child = (ChildRight(v) < heapSize) && g.compare(heap[ChildRight(v)], heap[ChildLeft(v)]) ? 
				    ChildRight(v) : ChildLeft(v);
            		if ( g.compare(x, heap[child]) ) break;
			else {
				heap[v] = heap[child];
				pos[heap[v]] = v;
				v = child;
			}
        	}
        	heap[v] = x;
		pos[x] = v;
    	}

	void initialize( GreaterActivity ga ) {
		g = ga;
	}

	int empty() { 
		if ( heapSize == 0 ) return 1;
       		else return 0;
	}

	int inHeap( int n ) { 
		if ( (n < posSize) && (pos[n] >= 0) ) return 1;
       		else return 0;	
	}

	void update( int x ) { up(pos[x]); }

    	void insert( int x ) {
		if ( posSize < x + 1 ) {
			for ( int i = posSize; i < x + 1; i ++ ) pos[i] = -1;
			posSize = x + 1;
		}
		pos[x] = heapSize;
		heap[heapSize] = x;
		heapSize++;
        	up(pos[x]); 
    	}

    	int pop() {
        	int x = heap[0];
        	heap[0] = heap[heapSize-1];
		pos[heap[0]] = 0;
		pos[x] = -1;
		heap[heapSize-1] = -1;
		heapSize--;
        	if ( heapSize > 1 ) down(0);
        	return x; 
    	}
} Heap;


// Clauses
class Clause {
public:
	// Literal block distance based on Glucose
	// LBD = How many decision variable in a learnt clause
    	int lbd;
    	int literals[40]; // Literals in this clause
	int literalsSize = 0;
    	int& operator [] ( int index ) { return literals[index]; } // Overloading array operator
    	Clause(): lbd(0) {}
	Clause( int sz ): lbd(0) { 
		if ( sz == literalsSize ) {
			literalsSize = sz;
		} else if ( sz > literalsSize ) {
			for ( int i = literalsSize; i < sz; i ++ ) {
				literals[i] = -1;
			}
			literalsSize = sz;
		} else {
			for ( int i = literalsSize; i < sz; i -- ) {
				literals[i] = -1;
			}
			literalsSize = sz;
		}
	} // Initialize int lbd as 0 and set a clause size with int size
	~Clause() {}
};


// The watched literals data structure (lazy data structure)
class WL {
public:
    	int clauseIdx; // Clause index in clause database
    	int blocker; // To fast guess whether a clause is already satisfied. 
    	WL(): clauseIdx(0), blocker(0) {}
    	WL( int c, int b ): clauseIdx(c), blocker(b) {}
};


// Solver
class Solver {
public:
    	int learnt[40]; 				// The clause indices of the learnt clauses
	int learntSize = 0;
	int trail[100]; 				// Save the assigned literal sequence(phase saving)
	int trailSize = 0;
	int decVarInTrail[50];	 			// Save the decision variables' position in trail(phase saving)
	int decVarInTrailSize = 0;
	int reduceMap[65000]; 				// Data structure for reduce
	int reduceMapSize = 0;
    	
	int vars, clauses, origin_clauses, conflicts;   // The number of variables, clauses, and conflicts
	int decides, propagations;			// The number of decides and propagations
    	int restarts, rephases, reduces;                // Parameters for restart, rephase, and reduce
    	int rephase_inc, rephase_limit, reduce_limit;   // Parameters for rephase and reduce
    	int threshold;                                  // A threshold for updating the local_best phase
    	int propagated;                                 // The number of propagted literals in trail
    	int time_stamp;                                 // Parameter for conflict analyzation and LBD calculation   
   
    	int lbd_queue[50],                              // Circled queue saved the recent 50 LBDs
            lbd_queue_size,                             // The number of LBDs in this queue
            lbd_queue_pos;                              // The position to save the next LBD
    	double fast_lbd_sum, slow_lbd_sum;              // Sum of the global and recent 50 LBDs

    	int *value,                                     // The variable assignement (1:True; -1:False; 0:Undefine) 
            *reason,                                    // The index of the clause that implies the variable assignment
            *level,                                     // The decision level of a variable      
            *mark,                                      // Parameter for conflict analyzation
            *local_best,                                // A phase with a local deepest trail                     
            *saved;                                     // Phase saving

    	double *activity;                               // The variables' score for VSIDS
    	double var_inc;                                 // Parameter for VSIDS
    	Heap vsids;                                     // Heap to select variable
     
    	void alloc_memory();                                      // Allocate memory 
    	void assign( int literal, int level, int cref );          // Assigned a variable
    	int  propagate();                                         // BCP (Boolean Contraint Propagation)
    	void backtrack( int backtrack_level );                    // Backtracking
    	int  analyze( int cref, int &backtrack_level, int &lbd ); // Conflict analyzation
    	int  parse( char *filename );                             // Read CNF file
	char *read_whitespace( char *p );
	char *read_until_new_line( char *p );
	char *read_int( char *p, int *i );
	int  solve();                                             // Solving
    	int  decide();                                            // Pick desicion variable
    	int  add_clause( int c[], int size );                     // Add new clause to clause database
    	void bump_var( int var, double mult );                    // Update activity      
    	void restart();                                           // Do restart                                      
    	void reduce();                                            // Do reduce
    	void rephase();                                           // Do rephase
    	void printModel();                                        // Print model when the result is SAT
};
