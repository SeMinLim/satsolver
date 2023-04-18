#include <stdbool.h>

#define ChildLeft(x) (x << 1 | 1)
#define ChildRight(x) ((x + 1) << 1)
#define Parent(x) ((x - 1) >> 1)

#define Value(literal) (literal > 0 ? value[literal] : -value[-literal])
#define WatchedLiterals(id) (watched_literals[s->vars + id])
#define WatchedLiteralsSize(id) (watched_literals_size[s->vars + id])

// Benchmark 7
#define NumVars 99
#define NumClauses 264
#define MaxNumLits 3


// Heap data structure (max heap)
typedef struct Heap {
	const double *activity; // Pointer to activity database
	int heap[NumVars+1]; // Index of activity array
	int heapSize;
	int pos[NumVars+1]; // Actual position of heap
	int posSize;
} Heap;
void heap_initialize( Heap *h, const double *a );
bool heap_compare( Heap *h, int a, int b );
void heap_up( Heap *h, int v );
void heap_down( Heap *h, int v );
bool heap_empty( Heap *h ); 
bool heap_inHeap( Heap *h, int n ); 
void heap_update( Heap *h, int x );
void heap_insert( Heap *h, int x );
int heap_pop( Heap *h );

		
// Clause
typedef struct Clause {
	// Literal Block Distance based on Glucose
	// LBD = How many decision variable in a learnt clause
    	int lbd;
    	int literals[64]; // Literals in a clause
	int literalsSize;
} Clause;
void clause_init( Clause *c );
void clause_resize( Clause *c, int sz ); // Init LBD as 0 and resize literal array


// Watcher list
typedef struct WL {
	// Which clause a watched literal is included
    	int clauseIdx; // A index of a clause in ClauseDB
    	int blocker; // A flag for check whether a clause is already satisfied
} WL;
void wl_init( WL *w );
void wl_set( WL *w, int c, int b );


// Solver
typedef struct Solver {
	int origin_clauses;				// # of origin clauses
	int vars, clauses, conflicts;			// # of variables, clauses, and conflicts
	int decides, propagations;			// # of decides and propagations
	int backtracklevel, lbd;			// Parameters for backtracking
    	int restarts;					// Parameter for restart
	int reduces, reduce_limit;			// Parameters for reduce
	int rephases, rephase_inc, rephase_limit;	// Parameters for rephase
    	int threshold;					// A threshold to update local-best phase
    	int propagated;					// # of propagated literals in trail
	int time_stamp;					// Parameter for conflict analysis

    	int lbd_queue[50],  	// Circles queue saving the recent 50 LBDs
            lbd_queue_size, 	// The number of LBDs in this queue
            lbd_queue_pos;  	// The position to save the next LBD
    	int fast_lbd_sum,	// Sum of the global LBDs 
	    slow_lbd_sum;	// Sum of the recent 50 LBDs
} Solver;
void solver_init( Solver *s );
void solver_assign( Solver *s, int literal, int level, int cref );
int  solver_add_clause( Solver *s, int c[], int size );
int  solver_propagate( Solver *s );
int  solver_parse( Solver *s );
int  solver_decide( Solver *s );
void solver_update_score( Solver *s, int var, double coeff );
int  solver_analyze( Solver *s, int cref );
void solver_backtrack( Solver *s, int backtrack_level );
void solver_restart( Solver *s );
void solver_reduce( Solver *s );
void solver_rephase( Solver *s );
int  solver_solve( Solver *s );


// Etc
unsigned int rand_generator();
int abs_value( int n );
void *mem_cpy( Clause *dest, Clause *src );
