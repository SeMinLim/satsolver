#include "solver.h"


// Global variables
int learnt[512];
int learntSize = 0;
int reduceMap[64*1024];
int reduceMapSize = 0;
Clause clauseDB[512*1024];
int clauseDBSize = 0;
WL watchedPointers[NumVars*2+1][32*1024];
int watchedPointersSize[NumVars*2+1] = {0,};


// Heap data structure
void heap_initialize( Heap *h, const double *a ) {
	h->activity = a;
	h->heapSize = 0;
	h->posSize = 0;
}

int heap_compare( Heap *h, int a, int b ) {
	if ( h->activity[a] > h->activity[b] ) return 1; 
	else return 0;
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

int heap_empty( Heap *h ) { 
	if ( h->heapSize == 0 ) return 1;
	else return 0;
}

int heap_inHeap( Heap *h, int n ) { 
	if ( (n < h->posSize) && (h->pos[n] >= 0) ) return 1;
	else return 0;	
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


// Solver
char *Solver::read_whitespace( char *p ) {
        while ( (*p >= 9 && *p <= 13) || *p == 32 ) ++p;
        return p;
}

char *Solver::read_until_new_line( char *p ) {
        while ( *p != '\n' ) {
                if ( *p++ == '\0' ) exit(1);
        }
        return ++p;
}

char *Solver::read_int( char *p, int *i ) {
        int sym = 1;
        *i = 0;
        p = read_whitespace(p);
        if ( *p == '-' ) {
                sym = 0;
                ++p;
        }
        while ( *p >= '0' && *p <= '9' ) {
                if ( *p == '\0' ) return p;
                *i = *i * 10 + *p - '0';
                ++p;
        }
        if ( !sym ) *i = -(*i);
        return p;
}

int Solver::add_clause( int c[], int size ) {
	Clause cls;
	cls.resize(size);
    	clauseDB[clauseDBSize++] = cls;
	int id = clauseDBSize - 1;                                
    	for ( int i = 0; i < size; i++ ) clauseDB[id].literals[i] = c[i];
	WL wl;
	wl.set(id, c[1]);
    	WatchedPointers(-c[0])[WatchedPointersSize(-c[0])++] = wl;
	wl.set(id, c[0]);
	WatchedPointers(-c[1])[WatchedPointersSize(-c[1])++] = wl;
	return id;                                                      
}

int Solver::parse( char *filename ) {
    	FILE *f_data = fopen(filename, "r");  

    	fseek(f_data, 0, SEEK_END);
    	size_t file_len = ftell(f_data);

	fseek(f_data, 0, SEEK_SET);
	char *data = new char[file_len + 1];
	char *p = data;
	fread(data, sizeof(char), file_len, f_data);
	fclose(f_data);                                             
	data[file_len] = '\0';

	int buffer[3];
	int bufferSize = 0;
	while ( *p != '\0' ) {
        	p = read_whitespace(p);
        	if ( *p == '\0' ) break;
        	if ( *p == 'c' ) p = read_until_new_line(p);
       	 	else if ( *p == 'p' ) { 
            		if ( (*(p + 1) == ' ') && (*(p + 2) == 'c') && 
			     (*(p + 3) == 'n') && (*(p + 4) == 'f') ) {
                		p += 5; 
				p = read_int(p, &vars); 
				p = read_int(p, &clauses);
                		alloc_memory();
            		} 
            		else printf("PARSE ERROR(Unexpected Char)!\n"), exit(2);
        	}
        	else {                                                                             
            		int32_t dimacs_lit;
            		p = read_int(p, &dimacs_lit);
            		if ( dimacs_lit != 0 ) { 
				if ( *p == '\0' ) {
                			printf("c PARSE ERROR(Unexpected EOF)!\n");
					exit(1);
				}
				else buffer[bufferSize++] = dimacs_lit;
			}
			else {                                                       
                		if ( bufferSize == 0 ) return 20;
				else if ( bufferSize == 1 ) {
					if ( Value(buffer[0]) == -1 ) return 20;
					else if ( !Value(buffer[0]) ) assign(buffer[0], 0, -1);
				}
                		else add_clause(buffer, bufferSize);
				bufferSize = 0;
            		}
        	}
    	}
    	origin_clauses = clauseDBSize;
    	return ( propagate() == -1 ? 0 : 20 );             
}

void Solver::alloc_memory() {
	conflicts = decides = propagations = restarts = rephases = reduces = 0;
	threshold = propagated = time_stamp = 0;
    	fast_lbd_sum = lbd_queue_size = lbd_queue_pos = slow_lbd_sum = 0;
    	var_inc = 1, rephase_inc = 1e5, rephase_limit = 1e5, reduce_limit = 8192;

    	// Initialization
	heap_initialize(&vsids, activity);
    	for (int i = 1; i <= vars; i++) {
        	value[i] = reason[i] = level[i] = mark[i] = local_best[i] = activity[i] = saved[i] = 0;
		heap_insert(&vsids, i);
    	}
}

void Solver::assign( int literal, int l, int cref ) {
    	int var = abs(literal);
    	value[var]  = literal > 0 ? 1 : -1;
    	level[var]  = l;
	reason[var] = cref;
	trail[trailSize++] = literal;
}

int Solver::decide() {      
    	int next = -1;
	// Pick up a variable based on VSIDS
	while ( next == -1 || Value(next) != 0 ) {
        	if (heap_empty(&vsids)) return 10;
        	else next = heap_pop(&vsids);
    	}
	decVarInTrail[decVarInTrailSize++] = trailSize;
    	// If there's saved one(polarity), use that
	if ( saved[next] ) next *= saved[next];

    	assign(next, decVarInTrailSize, -1);
    	decides++;
	return 0;
}

int Solver::propagate() {
	// This propagate style is fully based on MiniSAT
    	while ( propagated < trailSize ) { 
        	int p = trail[propagated++]; // 'p' is enqueued fact to propagate
		int size = WatchedPointersSize(p);
        	int i, j;                     
        	for ( i = j = 0; i < size;  ) {
			// Try to avoid inspecting the clause
            		int blocker = WatchedPointers(p)[i].blocker;                       
			if ( Value(blocker) == 1 ) {                
                		WatchedPointers(p)[j++] = WatchedPointers(p)[i++]; 
				continue;
            		}
            		// Make sure the false literal is 'c[1]'
			int cref = WatchedPointers(p)[i].clauseIdx;
			int falseLiteral = -p;
            		Clause& c = clauseDB[cref];              
            		if ( c.literals[0] == falseLiteral ) {
				c.literals[0] = c.literals[1];
				c.literals[1] = falseLiteral;
			}
			i++;
			// If 0th watch pointer is true, then clause is already satisfied
			int firstWP = c.literals[0];
			WL w;
		       	w.set(cref, firstWP);
			if ( Value(firstWP) == 1 ) {                   
                		WatchedPointers(p)[j++] = w; 
				continue;
            		}
			// Look for new watch pointer
			int k;
			int sz = c.literalsSize;
            		for ( k = 2; (k < sz) && (Value(c.literals[k]) == -1); k++ ); 
			if ( k < sz ) {                           
                		c.literals[1] = c.literals[k];
				c.literals[k] = falseLiteral;
				WatchedPointers(-c.literals[1])[WatchedPointersSize(-c.literals[1])++] = w;
			}         
			else { // Did not find new watch, clause is unit under assignment
				WatchedPointers(p)[j++] = w;
				// Conflict!
                		if ( Value(firstWP) == -1 ) { 
                    			while ( i < size ) WatchedPointers(p)[j++] = WatchedPointers(p)[i++];
					// Shrink
					WatchedPointersSize(p) = j;
                    			return cref;
                		}
				// Not conflict!
                		else {
					assign(firstWP, level[abs(p)], cref);
					propagations++;
				}
			}
            	}
		// Shrink
		WatchedPointersSize(p) = j;
    	}
    	return -1;                                       
}

void Solver::bump_var( int var, double coeff ) {
	// Update score and prevent float overflow
    	if ( (activity[var] += var_inc * coeff) > 1e100 ) {
        	for ( int i = 1; i <= vars; i++ ) activity[i] *= 1e-100;
        	var_inc *= 1e-100;
	}
	// Update Heap
    	if ( heap_inHeap(&vsids, var) ) heap_update(&vsids, var);
}

int Solver::analyze( int conflict, int &backtrackLevel, int &lbd ) {
    	++time_stamp;
    	learntSize = 0;
    	Clause &c = clauseDB[conflict]; 
	int conflictLevel = level[abs(c.literals[0])];

    	if ( conflictLevel == 0 ) return 20; // UNSAT
	else {
		learnt[learntSize++] = 0; // Leave a place to save the first UIP
		int bump[128];
		int bumpSize = 0;
		int should_visit_ct = 0; // The number of literals that have not visited in the conflict level of the implication graph
		int resolve_lit = 0; // The literal to do resolution
		int index = trailSize - 1;
		// First UIP learning method
		do {
			Clause &c = clauseDB[conflict];
			// Mark the literals
			for ( int i = (resolve_lit == 0 ? 0 : 1); i < c.literalsSize; i++ ) {
				int var = abs(c.literals[i]);
				if ( mark[var] != time_stamp && level[var] > 0 ) {
					bump_var(var, 0.5);
					bump[bumpSize++] = var;
					mark[var] = time_stamp;
					if ( level[var] >= conflictLevel ) should_visit_ct++;
					else learnt[learntSize++] = c.literals[i];
				}
			}
			// Find the last marked literal in the trail to do resolution
			do {
				while ( mark[abs(trail[index--])] != time_stamp );
				resolve_lit = trail[index + 1];
			} while ( level[abs(resolve_lit)] < conflictLevel );
			
			conflict = reason[abs(resolve_lit)];
			mark[abs(resolve_lit)] = 0;
			should_visit_ct--;
		} while ( should_visit_ct > 0 );

		learnt[0] = -resolve_lit;
		++time_stamp;
		lbd = 0;
		
		// Calculate LBD
		for ( int i = 0; i < learntSize; i++ ) {
			int l = level[abs(learnt[i])];
			if ( l && mark[l] != time_stamp ) {
				mark[l] = time_stamp;
				++lbd;
			}
		}

		if ( lbd_queue_size < 50 ) lbd_queue_size++;
		else fast_lbd_sum -= lbd_queue[lbd_queue_pos];
		
		fast_lbd_sum += lbd; // Sum of the recent 50 LBDs
		lbd_queue[lbd_queue_pos++] = lbd;
		
		if ( lbd_queue_pos == 50 ) lbd_queue_pos = 0;
		slow_lbd_sum += (lbd > 50 ? 50 : lbd); // Sum of the global LBDs
			
		if ( learntSize == 1 ) backtrackLevel = 0;
		else {
			int max_id = 1;
			for ( int i = 2; i < learntSize; i++ ) {
				if ( level[abs(learnt[i])] > level[abs(learnt[max_id])] ) max_id = i;
			}
			int p = learnt[max_id];
			learnt[max_id] = learnt[1];
			learnt[1] = p;
			backtrackLevel = level[abs(p)];
		}

		for ( int i = 0; i < bumpSize; i++ ) {   
			if ( level[bump[i]] >= backtrackLevel - 1 ) bump_var(bump[i], 1);
		}
	}
    	return 0;
}

void Solver::backtrack( int backtrackLevel ) {
    	if ( decVarInTrailSize <= backtrackLevel ) return;
	else {
		for ( int i = trailSize - 1; i >= decVarInTrail[backtrackLevel]; i-- ) {
			int v = abs(trail[i]);
			value[v] = 0;
			saved[v] = trail[i] > 0 ? 1 : -1; // Phase saving
			if ( !heap_inHeap(&vsids, v) ) heap_insert(&vsids, v); // Update heap
		}
		propagated = decVarInTrail[backtrackLevel];
		// Resize 'trail'
		for ( int i = trailSize; i < propagated; i -- ) trail[i] = -1;
		trailSize = propagated;
		// Resize 'decVarInTrail'
		for ( int i = decVarInTrailSize; i < backtrackLevel; i -- ) decVarInTrail[i] = -1;
		decVarInTrailSize = backtrackLevel;
	}
}

void Solver::restart() {
    	fast_lbd_sum = lbd_queue_size = lbd_queue_pos = 0;
    	backtrack(decVarInTrailSize);
	restarts++;
}

void Solver::rephase() {
	if ( (rephases / 2) == 1 ) for ( int i = 1; i <= vars; i++ ) saved[i] = local_best[i];
	else for ( int i = 1; i <= vars; i++ ) saved[i] = -local_best[i];
	backtrack(decVarInTrailSize);
	rephase_inc *= 2;
	rephase_limit = conflicts + rephase_inc;
	rephases++;
}

void Solver::reduce() {
    	backtrack(0);
    	reduces = 0;
	reduce_limit += 512;
    	int new_size = origin_clauses;
	int old_size = clauseDBSize;
    	// Resize 'reduceMap'
	if ( old_size > reduceMapSize ) {
		for ( int i = reduceMapSize; i < old_size; i ++ ) reduceMap[i] = -1;
		reduceMapSize = old_size;
	} else {
		for ( int i = reduceMapSize; i < old_size; i -- ) reduceMap[i] = -1;
		reduceMapSize = old_size;
	}
	// Random delete 50% bad clauses (LBD>=5) 
	// Reducing based on Literal Block Distances
    	for ( int i = origin_clauses; i < old_size; i++ ) { 
        	if ( clauseDB[i].lbd >= 5 && rand() % 2 == 0 ) reduceMap[i] = -1; // remove clause
        	else {
            		if ( new_size != i ) clauseDB[new_size] = clauseDB[i];
            		reduceMap[i] = new_size++;
        	}
    	}
    	//clauseDB.resize(new_size, Clause(0));
	if ( new_size > clauseDBSize ) {
		for ( int i = clauseDBSize; i < new_size; i ++ ) {
			Clause cls_1;
			cls_1.resize(0);
			clauseDB[i] = cls_1;
		}
	} else {
		for ( int i = clauseDBSize; i < new_size; i -- ) {
			Clause cls_2;
			cls_2.resize(0);
			clauseDB[i] = cls_2;
		}
	}
	clauseDBSize = new_size;

	for ( int v = -vars; v <= vars; v++ ) { // Update Watched Pointers
        	if ( v == 0 ) continue;
        	int old_sz = WatchedPointersSize(v);
		int new_sz = 0;

        	for ( int i = 0; i < old_sz; i++ ) {
            		int old_idx = WatchedPointers(v)[i].clauseIdx;
            		int new_idx = old_idx < origin_clauses ? old_idx : reduceMap[old_idx];
            		if ( new_idx != -1 ) {
                		WatchedPointers(v)[i].clauseIdx = new_idx;
                		if (new_sz != i) WatchedPointers(v)[new_sz] = WatchedPointers(v)[i];
                		new_sz++;
            		}
        	}
        	//WatchedPointers(v).resize(new_sz);
		WatchedPointersSize(v) = new_sz;
    	}
}

int Solver::solve() {
    	int res = 0;
    	while (!res) {
		int cref = propagate(); // Boolean Constraint Propagation (BCP)
		// Find a conflict
		if ( cref != -1 ) {
			int backtrackLevel = 0; 
			int lbd = 0;
			res = analyze(cref, backtrackLevel, lbd); // Conflict analyze
			
			if ( res == 20 ) break; // Find a conflict in 0 decision level
			else {
				backtrack(backtrackLevel); // backtracking
				
				if ( learntSize == 1 ) assign(learnt[0], 0, -1); // Learnt a unit clause
				else {
					int cref = add_clause(learnt, learntSize); // Add a clause to data base.
					clauseDB[cref].lbd = lbd;
					assign(learnt[0], backtrackLevel, cref); // The learnt clause implies the assignment of the UIP variable.
				}
				var_inc *= (1 / 0.8); // var_decay for locality
				++conflicts, ++reduces;
				
				// Update the local-best phase
				if ( trailSize > threshold ) {
					threshold = trailSize;
					for ( int i = 1; i < vars + 1; i++ ) local_best[i] = value[i];
				}
			}
		}
		else if ( reduces >= reduce_limit ) reduce();
		else if ( (lbd_queue_size == 50) && (0.8*fast_lbd_sum/lbd_queue_size > slow_lbd_sum/conflicts) ) restart();
		else if ( conflicts >= rephase_limit ) rephase();
		else res = decide();
	}
	return res;
}

