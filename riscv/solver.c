#include "solver.h"


// Benchmark
const int benchmark[264][3] = {
	{14, 1, -74},{-38, -58, 46},{48, -67, 49},{2, -20, -85},
	{81, 64, -46},{72, 80, 12},{53, 78, 11},{-97, -43, 33},
	{89, -66, -53},{-92, -55, -27},{5, 87, -81},{48, 93, 8},
	{77, 5, 24},{17, -71, 50},{40, 33, 69},{70, 86, 34},
	{90, -72, 42},{-54, -15, 43},{-38, -91, 27},{-11, -78, 53},
	{23, -74, 95},{-79, 2, 75},{-35, -59, 80},{-46, -64, -81},
	{-95, 57, 28},{-78, -93, 75},{-19, -90, -4},{-9, 12, -89},
	{67, -48, 49},{31, -88, 18},{21, -94, 55},{42, 63, 22},
	{35, 80, 59},{13, -11, -76},{46, 38, 58},{45, -37, 23},
	{52, -82, 7},{-66, 14, 34},{-8, -96, -61},{-21, -55, -94},
	{-10, 47, -24},{-12, 9, -89},{-30, 73, 50},{-14, -1, -74},
	{81, -64, 46},{35, 25, 6},{-32, -16, 28},{79, 2, -75},
	{-64, 99, -71},{-90, -42, -72},{87, 19, 52},{51, -68, 37},
	{-28, -32, 16},{67, -49, 48},{78, -11, -53},{3, -65, -99},
	{-36, -22, -31},{57, 85, -92},{-55, 92, 27},{56, 49, -65},
	{51, 97, -6},{91, -38, -27},{16, 39, 9},{-23, 45, 37},
	{6, 51, -97},{43, 54, 15},{89, -9, -12},{25, 98, 84},
	{-63, -67, 41},{-87, 81, 5},{10, 47, 24},{63, -41, -67},
	{-13, -94, 82},{78, 75, 93},{33, -69, -40},{6, 97, -51},
	{60, -26, 44},{-51, -6, -97},{60, 62, 79},{-16, -9, 39},
	{40, 47, -41},{30, 50, -73},{76, 11, 13},{-65, 99, -3},
	{-79, -75, -2},{-18, 31, 88},{8, 61, -96},{-2, 20, -85},
	{-5, 81, 87},{21, -62, 36},{82, 94, 13},{77, 3, 73},
	{-53, 11, -78},{91, -20, 98},{-23, 74, 95},{86, -44, -59},
	{16, 32, 28},{-4, -76, -10},{-60, -44, -26},{52, -7, 82},
	{67, -41, -63},{-49, -48, -67},{-87, -81, -5},{-34, 66, 14},
	{-45, 37, 23},{-14, 74, 1},{64, -81, 46},{84, -98, -25},
	{-60, 44, 26},{-19, 52, -87},{-13, 76, -11},{72, -42, 90},
	{-1, 14, 74},{-14, 34, 66},{-58, 29, 88},{30, 54, 84},
	{-8, 96, 61},{3, 99, 65},{-97, -33, 43},{-52, -7, -82},
	{-17, 50, 71},{-34, 70, -86},{-85, 92, 57},{-26, -61, -32},
	{-99, 65, -3},{-95, -57, -28},{-37, -68, -51},{-13, 94, -82},
	{-35, 59, -80},{82, -52, 7},{20, -98, 91},{-34, 86, -70},
	{12, -72, -80},{1, -96, 15},{94, -55, 21},{53, -89, -66},
	{-84, -98, 25},{6, -25, -35},{58, -38, -46},{27, -92, 55},
	{-29, -39, 68},{29, -39, -68},{-20, -98, -91},{24, -77, -5},
	{-16, 9, -39},{-7, -83, -18},{87, -52, -19},{62, -60, -79},
	{61, -32, 26},{68, 39, 29},{19, -52, -87},{26, -61, 32},
	{-47, -24, 10},{18, 83, -7},{51, -37, 68},{-72, 80, -12},
	{-88, -31, -18},{12, 89, 9},{69, -33, -40},{-80, 72, -12},
	{31, -22, 36},{-43, -33, 97},{-13, -76, 11},{83, 56, 70},
	{66, 53, 89},{20, 85, 2},{-88, 29, 58},{79, -2, 75},
	{-29, 88, 58},{7, -83, 18},{-6, 25, -35},{43, 33, 97},
	{-68, -29, 39},{59, -44, -86},{64, 99, 71},{-34, -14, -66},
	{1, 96, -15},{-58, 38, -46},{-30, 54, -84},{65, 56, -49},
	{83, -56, -70},{-36, 21, 62},{69, -17, 45},{-5, -24, 77},
	{-91, 38, -27},{40, -33, -69},{22, 31, -36},{-71, -99, 64},
	{-42, -22, 63},{69, 17, -45},{37, -51, 68},{-47, -41, -40},
	{-30, -73, -50},{44, 86, 59},{-71, -50, -17},{-57, -92, -85},
	{83, -18, 7},{96, 15, -1},{34, -70, -86},{-31, 18, 88},
	{-45, -69, -17},{-21, -36, -62},{-70, -83, 56},{-96, -1, -15},
	{4, 19, -90},{-10, 76, 4},{-53, -89, 66},{-93, 8, -48},
	{3, -73, -77},{71, -64, -99},{98, -91, 20},{-37, -23, -45},
	{49, 65, -56},{-95, -74, -23},{73, -3, -77},{-84, 98, -25},
	{-69, 45, 17},{-86, -59, 44},{10, 76, -4},{-63, -22, 42},
	{-19, 4, 90},{-25, 35, -6},{-47, -10, 24},{-42, 22, -63},
	{-27, 55, 92},{-62, 79, -60},{-30, -54, 84},{-40, 47, 41},
	{19, 90, -4},{-61, 96, 8},{-48, 93, -8},{92, 85, -57},
	{-24, -77, 5},{-21, 36, 62},{-59, 35, -80},{57, -28, 95},
	{-75, -93, 78},{22, 36, -31},{41, 40, -47},{13, -94, -82},
	{55, 94, -21},{15, -43, -54},{-29, -58, -88},{60, 26, -44},
	{-54, 30, -84},{-65, -49, -56},{-78, 93, -75},{-15, -43, 54},
	{-83, 70, -56},{-50, 30, 73},{38, 27, 91},{-2, 85, -20},
	{67, 41, 63},{17, 71, -50},{60, -62, -79},{-95, 23, 74},
	{4, -76, 10},{16, -9, -39},{72, -90, 42},{77, -73, -3},
	{-28, 32, -16},{-57, 95, 28},{-26, 61, 32},{-8, 48, -93}
};


// Global variables
int learnt[512];				// The index of the learnt clauses
int learntSize = 0;				//
int trail[NumVars];				// Save the assigned literal sequence
int trailSize = 0;				// (phase saving)
int decVarInTrail[NumVars/2];			// Save the decision variables' position
int decVarInTrailSize = 0;			// in trail (phase saving)
int reduceMap[128*1024];			// Data structure for reduce
int reduceMapSize = 0;				//
Clause clauseDB[128*1024];			// Clause database
int clauseDBSize = 0;
WL watched_literals[NumVars*2+1][32*1024];	// A mapping from a literal to clauses
int watched_literals_size[NumVars*2+1] = {0,};	//
	
int8_t value[NumVars+1] = {0,};		// The variable assignment (1:True;-1:False;0:Undefine)
int8_t local_best[NumVars+1] = {0,};	// A phase with a local deepest trail
int8_t saved[NumVars+1] = {0,};		// Phase saving
int reason[NumVars+1] = {0,};		// The index of the clause that implies the variable assignment
int level[NumVars+1] = {0,};		// The decision level of a variable
int mark[NumVars+1] = {0,};		// Parameter for conflict analysis

uint64_t activity[NumVars+1] = {(uint64_t)0U,};	// The variables' score for VSIDS
Heap vsids;					// Heap to select variable


// Etc
// rand() in stdlib
uint16_t lfsr = 0xACE1u;
uint16_t rand_generator() {
	uint32_t bit  = ((lfsr >> 0) ^ (lfsr >> 2) ^ (lfsr >> 3) ^ (lfsr >> 5) ) & 1;
	lfsr =  (lfsr >> 1) | (bit << 15);
	uint16_t value = lfsr << 15;
	return value;
}

// abs() in stdlib
int abs_value( int n ){
	return (n + (n >> 31)) ^ (n >> 31);
}

// memcpy in stdlib
void *mem_cpy( Clause *dest, const Clause *src ) {
	dest->lbd = src->lbd;
	for ( int i = 0; i < 64; i ++ ) dest->literals[i] = src->literals[i];
	dest->literalsSize = src->literalsSize;
	return dest;
}


// Heap data structure (max heap)
void heap_initialize( Heap *h, const uint64_t *a ) {
	h->activity = a;
	h->heapSize = 0;
	h->posSize = 0;
}

bool heap_compare( Heap *h, int a, int b ) {
	return h->activity[a] > h->activity[b]; 
}

void heap_up( Heap *h, int v ) {
	int x = h->heap[v];
	int p = Parent(v);
        while ( v && heap_compare(h, x, h->heap[p]) ) {
       		h->heap[v] = h->heap[p];
		h->pos[h->heap[p]] = v;
		v = p; 
		p = Parent(p);
        }
        h->heap[v] = x;
	h->pos[x] = v;
}

void heap_down( Heap *h, int v ) {
	int x = h->heap[v];
        while ( ChildLeft(v) < h->heapSize ){
		int child = (ChildRight(v) < h->heapSize) && heap_compare(h, h->heap[ChildRight(v)], h->heap[ChildLeft(v)]) ? 
			    ChildRight(v) : ChildLeft(v);
		if ( heap_compare(h, x, h->heap[child]) ) break;
		else {
			h->heap[v] = h->heap[child];
			h->pos[h->heap[v]] = v;
			v = child;
		}
        }
        h->heap[v] = x;
	h->pos[x] = v;
}

bool heap_empty( Heap *h ) { 
	return h->heapSize == 0;
}

bool heap_inHeap( Heap *h, int n ) { 
	return n < h->posSize && h->pos[n] >= 0;
}

void heap_update( Heap *h, int x ) { heap_up(h, h->pos[x]); }

void heap_insert( Heap *h, int x ) {
	if ( h->posSize < x + 1 ) {
		for ( int i = h->posSize; i < x + 1; i ++ ) h->pos[i] = -1;
		h->posSize = x + 1;
	}
	h->pos[x] = h->heapSize;
	h->heap[h->heapSize] = x;
	h->heapSize++;
	heap_up(h, h->pos[x]); 
}

int heap_pop( Heap *h ) {
	int x = h->heap[0];
        h->heap[0] = h->heap[h->heapSize-1];
	h->pos[h->heap[0]] = 0;
	h->pos[x] = -1;
	h->heap[h->heapSize-1] = -1;
	h->heapSize--;
        if ( h->heapSize > 1 ) heap_down(h, 0);
        return x; 
}


// Clause
void clause_init( Clause *c ) { 
	c->lbd = 0;
	c->literalsSize = 0;
}
	
void clause_resize( Clause *c, int sz ) {
	c->lbd = 0;
	if ( sz == c->literalsSize ) {
		c->literalsSize = sz;
	} else if ( sz > c->literalsSize ) {
		for ( int i = c->literalsSize; i < sz; i ++ ) {
			c->literals[i] = -1;
		}
		c->literalsSize = sz;
	} else {
		for ( int i = c->literalsSize; i < sz; i -- ) {
			c->literals[i] = -1;
		}
		c->literalsSize = sz;
	}
}


// Watched Literals
void wl_init( WL *w ) {
	w->clauseIdx = 0;
	w->blocker = 0;
}

void wl_set( WL *w, int c, int b ) {
	w->clauseIdx = c;
	w->blocker = b;
}


// Solver
// Initialize the values
void solver_init( Solver *s ) {
	s->backtracklevel = s->lbd = 0;

	s->conflicts = s->decides = s->propagations = s->restarts = s->rephases = s->reduces = 0;
	s->threshold = s->propagated = s->time_stamp = 0;
    	s->fast_lbd_sum = s->lbd_queue_size = s->lbd_queue_pos = s->slow_lbd_sum = 0;
    	
	s->rephase_inc = 1e5, s->rephase_limit = 1e5, s->reduce_limit = 8192;

	heap_initialize(&vsids, activity);
    	for (int i = 1; i <= s->vars; i++) heap_insert(&vsids, i);
}

// Assign true value to a certain literal
void solver_assign( Solver *s, int literal, int l, int cref ) {
    	int var = abs_value(literal);
    	value[var]  = literal > 0 ? 1 : -1;
    	level[var]  = l;
	reason[var] = cref;
	trail[trailSize++] = literal;
}

// Add a clause to the database
int solver_add_clause( Solver *s, int c[], int size ) {
	Clause cls;
	clause_init(&cls);
	clause_resize(&cls, size);
	mem_cpy(&clauseDB[clauseDBSize], &cls);
	clauseDBSize++;

	int id = clauseDBSize - 1;                                
    	for ( int i = 0; i < size; i++ ) clauseDB[id].literals[i] = c[i];

	WL wl;
	wl_set(&wl, id, c[1]);
    	WatchedLiterals(-c[0])[WatchedLiteralsSize(-c[0])++] = wl;
	wl_set(&wl, id, c[0]);
	WatchedLiterals(-c[1])[WatchedLiteralsSize(-c[1])++] = wl;
	return id;                                                      
}

// BCP (Boolean Constraint Propagation)
int solver_propagate( Solver *s ) {
    	while ( s->propagated < trailSize ) { 
        	int p = trail[s->propagated++];
		
		WL *ws = WatchedLiterals(p);
		int num_clauses = WatchedLiteralsSize(p);
        	
		int new_num_clauses = 0;                     
        	for ( int i = 0; i < num_clauses;  ) {
            		int blocker = ws[i].blocker;                       
			if ( Value(blocker) == 1 ) {                
                		ws[new_num_clauses++] = ws[i++]; 
				continue;
            		}

			int cref = ws[i].clauseIdx;
			int falseLiteral = -p;
            		Clause *c = &clauseDB[cref];              
            		if ( c->literals[0] == falseLiteral ) {
				c->literals[0] = c->literals[1];
				c->literals[1] = falseLiteral;
			}
			
			i++;
			
			int firstWP = c->literals[0];
			WL w;
		       	wl_set(&w, cref, firstWP);
			if ( Value(firstWP) == 1 ) {                   
                		ws[new_num_clauses++] = w; 
				continue;
            		}
			
			int k;
			int sz = c->literalsSize;
            		for ( k = 2; (k < sz) && (Value(c->literals[k]) == -1); k++ ); 
			if ( k < sz ) {                           
                		c->literals[1] = c->literals[k];
				c->literals[k] = falseLiteral;
				WatchedLiterals(-c->literals[1])[WatchedLiteralsSize(-c->literals[1])++] = w;
			} else { 
				ws[new_num_clauses++] = w;
                		if ( Value(firstWP) == -1 ) { 
                    			while ( i < num_clauses ) ws[new_num_clauses++] = ws[i++];
					WatchedLiteralsSize(p) = new_num_clauses;
                    			return cref;
                		} else {
					solver_assign(s, firstWP, level[abs_value(p)], cref);
					s->propagations++;
				}
			}
            	}
		WatchedLiteralsSize(p) = new_num_clauses;
    	}
    	return -1;                                       
}

// Transfer data to ClauseDB
int solver_parse( Solver *s ) {
	int buffer[3];
	int bufferSize = 0;
	
	s->vars = NumVars;
	s->clauses = NumClauses;
	solver_init(s);
	
	for ( int i = 0; i < NumClauses; i ++) {
		for ( int j = 0; j < MaxNumLits; j ++ ) {
			buffer[j] = benchmark[i][j];
			bufferSize++;
		}
		solver_add_clause(s, buffer, bufferSize);
		bufferSize = 0;
	}

    	s->origin_clauses = clauseDBSize;

	return ( solver_propagate(s) == -1 ? 0 : 20 );
}

// Pick deicison variable based on VSIDS
int solver_decide( Solver *s ) {      
    	int next = -1;

	while ( next == -1 || Value(next) != 0 ) {
        	if (heap_empty(&vsids)) return 10;
        	else next = heap_pop(&vsids);
    	}
	decVarInTrail[decVarInTrailSize++] = trailSize;
	
	if ( saved[next] ) {
		if ( saved[next] < 0 ) next = -next;
		else if ( saved[next] > 0 ) next = next;
	}
	solver_assign(s, next, decVarInTrailSize, -1);
    	
	s->decides++;
	return 0;
}

// Update activity
void solver_update_score( Solver *s, int var, uint64_t amount ) {
    	// Update score and prevent overflow
	// Integer type bumping scheme
	if ( activity[var] + amount > (uint64_t)18446744073709551515U ) {
		for ( int i = 1; i <= s->vars; i ++ ) {
			activity[i] >>= 1;
		}
		activity[var] += amount>>1;
	} else activity[var] += amount;
    	if ( heap_inHeap(&vsids, var) ) heap_update(&vsids, var);
}

// Conflict analysis
int solver_analyze( Solver *s, int conflict ) {
    	++s->time_stamp;
    	learntSize = 0;
    	Clause *c = &clauseDB[conflict]; 
	int conflictLevel = level[abs_value(c->literals[0])];

    	if ( conflictLevel == 0 ) return 20;
	else {
		learnt[learntSize++] = 0;
		
		int bump[128];
		int bumpSize = 0;
		
		int should_visit_ct = 0; 
		int resolve_lit = 0;
		
		int index = trailSize - 1;
		do {
			Clause *c = &clauseDB[conflict];
			for ( int i = (resolve_lit == 0 ? 0 : 1); i < c->literalsSize; i++ ) {
				int var = abs_value(c->literals[i]);
				if ( mark[var] != s->time_stamp && level[var] > 0 ) {
					solver_update_score(s, var, (uint64_t)2U);
					bump[bumpSize++] = var;
					mark[var] = s->time_stamp;
					if ( level[var] >= conflictLevel ) should_visit_ct++;
					else learnt[learntSize++] = c->literals[i];
				}
			}
		
			do {
				while ( mark[abs_value(trail[index--])] != s->time_stamp );
				resolve_lit = trail[index + 1];
			} while ( level[abs_value(resolve_lit)] < conflictLevel );
			
			conflict = reason[abs_value(resolve_lit)];
			mark[abs_value(resolve_lit)] = 0;
			should_visit_ct--;
		} while ( should_visit_ct > 0 );

		learnt[0] = -resolve_lit;
		++s->time_stamp;
		s->lbd = 0;
		
		for ( int i = 0; i < learntSize; i++ ) {
			int l = level[abs_value(learnt[i])];
			if ( l && mark[l] != s->time_stamp ) {
				mark[l] = s->time_stamp;
				++s->lbd;
			}
		}

		if ( s->lbd_queue_size < 50 ) s->lbd_queue_size++;
		else s->fast_lbd_sum -= s->lbd_queue[s->lbd_queue_pos];
		
		s->fast_lbd_sum += s->lbd;
		s->lbd_queue[s->lbd_queue_pos++] = s->lbd;
		
		if ( s->lbd_queue_pos == 50 ) s->lbd_queue_pos = 0;
		s->slow_lbd_sum += (s->lbd > 50 ? 50 : s->lbd);
			
		if ( learntSize == 1 ) s->backtracklevel = 0;
		else {
			int max_id = 1;
			for ( int i = 2; i < learntSize; i++ ) {
				if ( level[abs_value(learnt[i])] > level[abs_value(learnt[max_id])] ) max_id = i;
			}
			int p = learnt[max_id];
			learnt[max_id] = learnt[1];
			learnt[1] = p;
			s->backtracklevel = level[abs_value(p)];
		}

		for ( int i = 0; i < bumpSize; i++ ) {   
			if ( level[bump[i]] >= s->backtracklevel - 1 ) solver_update_score(s, bump[i], (uint64_t)4U);
		}
	}
    	return 0;
}

// Backtracking
void solver_backtrack( Solver *s, int backtrackLevel ) {
    	if ( decVarInTrailSize <= backtrackLevel ) return;
	else {
		for ( int i = trailSize - 1; i >= decVarInTrail[backtrackLevel]; i-- ) {
			int v = abs_value(trail[i]);
			value[v] = 0;
			saved[v] = trail[i] > 0 ? 1 : -1;
			if ( !heap_inHeap(&vsids, v) ) heap_insert(&vsids, v);
		}
		s->propagated = decVarInTrail[backtrackLevel];
		
		trailSize = s->propagated;
		
		decVarInTrailSize = backtrackLevel;
	}
}

// Do restart
void solver_restart( Solver *s ) {
    	s->fast_lbd_sum = s->lbd_queue_size = s->lbd_queue_pos = 0;
    	solver_backtrack(s, decVarInTrailSize);
	s->restarts++;
}

// Do rephase
void solver_rephase( Solver *s ) {
	if ( (s->rephases >> 1) == 1 ) for ( int i = 1; i <= s->vars; i++ ) saved[i] = local_best[i];
	else for ( int i = 1; i <= s->vars; i++ ) saved[i] = -local_best[i];
	solver_backtrack(s, decVarInTrailSize);
	s->rephase_inc <<= 1;
	s->rephase_limit = s->conflicts + s->rephase_inc;
	s->rephases++;
}

// Do reduce
void solver_reduce( Solver *s ) {
    	solver_backtrack(s, 0);
    	
	s->reduces = 0;
	s->reduce_limit += 512;
    	
	int new_size = s->origin_clauses;
	int old_size = clauseDBSize;

	reduceMapSize = old_size;

    	for ( int i = s->origin_clauses; i < old_size; i++ ) { 
        	if ( clauseDB[i].lbd >= 5 && rand_generator() == 0 ) reduceMap[i] = -1;
        	else {
            		if ( new_size != i ) {
				mem_cpy(&clauseDB[new_size], &clauseDB[i]);
			}
            		reduceMap[i] = new_size++;
        	}
    	}

	for ( int i = clauseDBSize; i < new_size; i -- ) {
		Clause cls_2;
		clause_resize(&cls_2, 0);
		mem_cpy(&clauseDB[i], &cls_2);
	}	
	clauseDBSize = new_size;

	for ( int v = -s->vars; v <= s->vars; v++ ) {
        	if ( v == 0 ) continue;
        
		int old_sz = WatchedLiteralsSize(v);
		int new_sz = 0;

		WL *ws =  WatchedLiterals(v);
        	for ( int i = 0; i < old_sz; i++ ) {
            		int old_idx = ws[i].clauseIdx;
            		int new_idx = old_idx < s->origin_clauses ? old_idx : reduceMap[old_idx];
            		if ( new_idx != -1 ) {
                		ws[i].clauseIdx = new_idx;
                		if (new_sz != i) ws[new_sz] = ws[i];
                		new_sz++;
            		}
        	}
		WatchedLiteralsSize(v) = new_sz;
    	}
}

// Solver
int solver_solve( Solver *s ) {
    	int res = 0;
    	while (!res) {
		int cref = solver_propagate(s);

		if ( cref != -1 ) {
			s->backtracklevel = 0; 
			s->lbd = 0;
			res = solver_analyze(s, cref);
			
			if ( res == 20 ) break;
			else {
				solver_backtrack(s, s->backtracklevel);
				
				if ( learntSize == 1 ) solver_assign(s, learnt[0], 0, -1);
				else {
					int cref = solver_add_clause(s, learnt, learntSize);
					clauseDB[cref].lbd = s->lbd;
					solver_assign(s, learnt[0], s->backtracklevel, cref);
				}
				// var_decay for locality
				for ( int i = 1; i <= s->vars; i ++ ) if ( activity[i] != 0 ) activity[i] -= (uint64_t)1U;

				++s->conflicts, ++s->reduces;
				
				if ( trailSize > s->threshold ) {
					s->threshold = trailSize;
					for ( int i = 1; i < s->vars + 1; i++ ) local_best[i] = value[i];
				}
			}
		} else if ( s->reduces >= s->reduce_limit ) {
			solver_reduce(s);
		} else if ( s->lbd_queue_size == 50 && 
			    s->fast_lbd_sum/s->lbd_queue_size > s->slow_lbd_sum/s->conflicts ) {
			solver_restart(s);
		} else if ( s->conflicts >= s->rephase_limit ) {
			solver_rephase(s);
		} else res = solver_decide(s);
	}
	return res;
}

