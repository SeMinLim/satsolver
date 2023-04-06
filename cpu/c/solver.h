#include <stdio.h>
#include <stdlib.h>

#define ChildLeft(x) (x << 1 | 1)
#define ChildRight(x) ((x + 1) << 1)
#define Parent(x) ((x - 1) >> 1)

#define Value(literal) (literal > 0 ? s->value[literal] : -s->value[-literal])
#define WatchedPointers(id) (watchedPointers[s->vars + id])
#define WatchedPointersSize(id) (watchedPointersSize[s->vars + id])

#define NumVars 99
#define NumClauses 264
#define MaxNumLits 3


// Heap data structure
typedef struct Heap {
	const double *activity;
	int heap[NumVars+1];
	int heapSize;
	int pos[NumVars+1];
	int posSize;
} Heap;
void heap_initialize( Heap *h, const double *a );
int heap_compare( Heap *h, int a, int b );
void heap_up( Heap *h, int v );
void heap_down( Heap *h, int v );
int heap_empty( Heap *h ); 
int heap_inHeap( Heap *h, int n ); 
void heap_update( Heap *h, int x );
void heap_insert( Heap *h, int x );
int heap_pop( Heap *h );

		
// Clauses
typedef struct Clause {
    	int lbd;
    	int literals[64];
	int literalsSize;
} Clause;
void clause_init( Clause *c );
void clause_resize( Clause *c, int sz );


// The watched literals data structure (lazy data structure)
typedef struct WL {
    	int clauseIdx;
    	int blocker;
} WL;
void wl_init( WL *w );
void wl_set( WL *w, int c, int b );


// Solver
typedef struct Solver {
	int trail[NumVars];
	int trailSize;
	int decVarInTrail[NumVars/2];
	int decVarInTrailSize;

	int vars, clauses, origin_clauses, conflicts;	
	int decides, propagations;
	int backtracklevel, lbd;
    	int restarts, rephases, reduces;
    	int rephase_inc, rephase_limit, reduce_limit;
    	int threshold;
    	int propagated;
	int time_stamp;
   
    	int lbd_queue[50],
            lbd_queue_size,
            lbd_queue_pos;
    	double fast_lbd_sum, slow_lbd_sum;

	int value[NumVars+1];
	int reason[NumVars+1];
	int level[NumVars+1];
	int mark[NumVars+1];
	int local_best[NumVars+1];
	int saved[NumVars+1];

	double activity[NumVars+1];
    	double var_inc;
    	Heap vsids; 
} Solver;
void solver_alloc_memory( Solver *s );
void solver_assign( Solver *s, int literal, int level, int cref );
int  solver_propagate( Solver *s );
void solver_backtrack( Solver *s, int backtrack_level );
int  solver_analyze( Solver *s, int cref );
int  solver_parse( Solver *s );
char *read_whitespace( char *p );
char *read_until_new_line( char *p );
char *read_int( char *p, int *i );
int  solver_solve( Solver *s );
int  solver_decide( Solver *s );
int  solver_add_clause( Solver *s, int c[], int size );
void solver_bump_var( Solver *s, int var, double mult );
void solver_restart( Solver *s );
void solver_reduce( Solver *s );
void solver_rephase( Solver *s );
