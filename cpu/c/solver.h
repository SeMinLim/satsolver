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
    	int lbd;
    	int literals[64];
	int literalsSize = 0;
    	int& operator [] ( int index ) { return literals[index]; }
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
	}
	~Clause() {}
};


// The watched literals data structure (lazy data structure)
class WL {
public:
    	int clauseIdx;
    	int blocker;
    	WL(): clauseIdx(0), blocker(0) {}
    	WL( int c, int b ): clauseIdx(c), blocker(b) {}
};


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

    	int *value,
            *reason,
            *level,
            *mark,
            *local_best,
            *saved;

    	double *activity;
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
