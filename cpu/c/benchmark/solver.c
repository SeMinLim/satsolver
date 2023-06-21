#include "solver.h"


// Benchmark
const int benchmark[218][3] = { 
	{-48, -49, 21},{41, 37, -42},{-25, -22, 10},{-25, 26, -8},
	{-20, 27, -4},{2, 9, 36},{-44, 47, 10},{-29, 28, -33},
	{6, -36, 8},{3, -49, -38},{-43, -27, 8},{2, -42, 9},
	{23, -5, -39},{-37, 49, 36},{-5, 16, -13},{-40, 14, 47},
	{44, -4, -18},{-28, -35, 19},{19, -26, 4},{18, 4, 30},
	{-50, -38, -31},{6, -24, 46},{3, -6, -10},{-40, 7, 12},
	{-11, -37, -4},{-8, -45, -16},{-31, 47, -19},{-2, -23, -11},
	{-34, 46, 9},{-37, -1, -13},{34, -42, 6},{-3, -29, 20},
	{-31, 21, 19},{-37, 48, 3},{4, -10, 18},{-6, -21, 12},
	{46, 17, -33},{-5, -35, 18},{-32, 46, -28},{30, -13, 45},
	{-37, 33, 12},{-25, 37, 7},{-3, -36, -23},{29, -30, 25},
	{-34, -29, 31},{-29, 32, 40},{14, -36, -47},{-26, -37, 33},
	{46, -13, -20},{-1, 24, -19},{-9, 49, 15},{-12, -45, 21},
	{8, -16, 50},{-1, 12, 21},{-30, 27, 10},{46, -17, -16},
	{-49, -13, 9},{-43, -30, 22},{-4, 6, 43},{13, -42, -28},
	{40, -34, -47},{-32, 37, -1},{-6, 9, -27},{-1, 10, -41},
	{46, 13, -20},{8, 34, -46},{-44, 8, 34},{39, 26, -23},
	{15, 24, 36},{47, 3, 20},{15, -27, 4},{44, 32, -28},
	{4, 28, -13},{-34, -25, -39},{33, 1, -3},{-25, 26, 19},
	{-30, -47, 16},{10, -1, 39},{35, -14, -11},{49, 45, 47},
	{-31, -40, 38},{50, -18, -4},{34, 10, -41},{45, 6, -43},
	{-50, 20, -32},{-8, 2, -38},{2, 1, -44},{-45, 49, 3},
	{15, 16, -28},{-45, 20, -18},{24, -39, -11},{19, 40, 30},
	{18, -4, -25},{17, 41, 25},{-18, -31, -41},{-18, 20, 29},
	{-7, 47, 44},{-47, -6, -40},{-16, -21, -12},{8, 12, 26},
	{-19, 28, -22},{11, 21, -26},{20, -47, 33},{46, 12, 34},
	{42, 8, -27},{-22, 28, -13},{-22, -36, -49},{25, 9, -46},
	{19, 13, 2},{-17, 11, -44},{-50, 47, 41},{15, -7, 17},
	{-37, 5, 11},{12, 34, 3},{31, -13, -16},{-3, -41, 32},
	{-12, 49, 9},{43, -47, -9},{4, 2, 13},{-11, 46, 29},
	{-11, 22, -10},{45, 26, 42},{8, -46, -49},{47, 2, 40},
	{25, -19, 50},{-22, -5, 41},{25, -41, -30},{19, 2, -28},
	{-6, 15, -27},{-50, 48, -2},{18, 23, -1},{-12, -29, -37},
	{24, -14, -32},{27, -36, -14},{-37, -30, -35},{-35, 46, 8},
	{32, 16, -40},{-49, -31, -48},{-30, 8, -15},{31, -47, 25},
	{1, 20, 50},{19, -36, -8},{-2, -6, 41},{13, 45, -34},
	{37, -3, 4},{3, 7, -11},{-46, 28, -37},{-1, -9, -19},
	{45, -42, -21},{-38, -49, 4},{42, 43, -17},{31, -50, 4},
	{46, -49, 16},{2, -27, -3},{-23, 15, -48},{37, -45, -26},
	{-30, -27, -29},{40, 25, 5},{-6, 33, 18},{-1, -24, -13},
	{13, 22, -17},{-45, 30, 29},{-28, -31, -23},{37, -20, 5},
	{-15, -3, -34},{4, 2, 21},{14, -32, 45},{30, 18, 35},
	{48, 25, -27},{26, -7, 31},{-21, 38, 12},{7, 21, -46},
	{-14, 3, 1},{-20, 11, -30},{44, 48, -41},{-35, 14, 47},
	{-32, 41, -6},{-47, -17, -39},{14, 42, 24},{-4, -11, 48},
	{13, -11, 38},{20, 18, 49},{14, -47, 29},{46, 5, 50},{27, -19, -24},
	{48, -9, -35},{38, 12, 40},{-26, -9, -30},{-23, 4, -38},
	{-19, -38, 24},{40, 43, 31},{15, 48, -36},{6, 35, 36},
	{-35, -9, 11},{-46, -15, 34},{-27, 41, 36},{15, 49, -47},
	{26, 44, 50},{-34, 40, -14},{-5, -12, 6},{-44, 25, -45},
	{40, 44, -37},{-34, -45, -30},{28, 39, -48},{27, 16, -44},
	{-8, 39, -42},{27, 3, -30},{43, -44, -27},{-8, 15, -5},
	{-13, -36, -37},{35, 39, 26},{48, 5, 42},{-4, -9, 49},
	{-33, -11, 39},{32, -10, -44},{-3, -22, -33},{-39, 42, 40},
	{3, -5, 14} 
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
int watched_literals_size[NumVars*2+1];		//
	
int8_t value[NumVars+1];		// The variable assignment (1:True;-1:False;0:Undefine)
int8_t local_best[NumVars+1];		// A phase with a local deepest trail
int8_t saved[NumVars+1];		// Phase saving
int reason[NumVars+1];			// The index of the clause that implies the variable assignment
int level[NumVars+1];			// The decision level of a variable
int mark[NumVars+1];			// Parameter for conflict analysis

uint64_t activity[NumVars+1];		// The variables' score for VSIDS
Heap vsids;				// Heap to select variable


// Etc
// rand() in stdlib (AND, Shift, XOR)
uint32_t lfsr32 = 0xACE8F;
uint32_t lfsr31 = 0x23456789;
uint32_t shift_lfsr( uint32_t *lfsr, uint32_t polynomial_mask ) {
	uint32_t feedback = *lfsr & 1;
	*lfsr >>= 1;
	if ( feedback == 1 ) *lfsr ^= polynomial_mask;
	return *lfsr;
}
uint32_t rand_generator() {
	shift_lfsr(&lfsr32, POLY_MASK_32);
	uint32_t tmp_1 = (shift_lfsr(&lfsr32, POLY_MASK_32) ^ shift_lfsr(&lfsr31, POLY_MASK_31)) & 0xffff;
	uint32_t tmp_2 = tmp_1 << 31;
	uint32_t value = tmp_2 >> 31;
	return value;
}

// abs() in stdlib (Shift, XOR)
int abs_value( int n ){
	return (n + (n >> 31)) ^ (n >> 31);
}

// memcpy in stdlib
void *mem_cpy( Clause *dest, const Clause *src ) {
	dest->lbd = src->lbd;
	for ( int i = 0; i < MaxNumLits; i ++ ) dest->literals[i] = src->literals[i];
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
	for ( int i = 0; i < MaxNumLits; i ++ ) {
		c->literals[i] = 0;
	}
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
		clause_init(&cls_2);
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

