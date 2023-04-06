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
	int literalsSize = 0;
	void init() { lbd = 0; }
	void resize( int sz ) {
		lbd = 0;
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
	}
} Clause;


// The watched literals data structure (lazy data structure)
typedef struct WL {
    	int clauseIdx;
    	int blocker;
	void init() {
		clauseIdx = 0;
		blocker = 0;
	}
    	void set( int c, int b ) {
		clauseIdx = c;
		blocker = b;
	}
} WL;


// Solver
class Solver {
public:
	int trail[NumVars];
	int trailSize = 0;
	int decVarInTrail[NumVars/2];
	int decVarInTrailSize = 0;

	int vars, clauses, origin_clauses, conflicts;	
	int decides, propagations;
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
     
    	void alloc_memory();
    	void assign( int literal, int level, int cref );
    	int  propagate();
    	void backtrack( int backtrack_level );
    	int  analyze( int cref, int &backtrack_level, int &lbd );
    	int  parse( char *filename );
	char *read_whitespace( char *p );
	char *read_until_new_line( char *p );
	char *read_int( char *p, int *i );
	int  solve();
    	int  decide();
    	int  add_clause( int c[], int size );
    	void bump_var( int var, double mult );
    	void restart();
    	void reduce();
    	void rephase();
};
