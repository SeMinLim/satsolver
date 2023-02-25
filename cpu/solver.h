#include <vector>
#include <fstream>

#define ChildLeft(x) (x << 1 | 1)
#define ChildRight(x) ((x + 1) << 1)
#define Parent(x) ((x - 1) >> 1)

// Heap Data Structure
template<class Compare>
class Heap {
    	Compare lt;
    	std::vector<int> heap;
    	std::vector<int> pos;
    
    	void up(int v) {
        	int x = heap[v];
		int p = Parent(v);
        	while (v && lt(x, heap[p])) {
       			heap[v] = heap[p];
			pos[heap[p]] = v;
            		v = p; 
			p = Parent(p);
        	}
        	heap[v] = x;
		pos[x] = v;
    	}

    	void down(int v) {
        	int x = heap[v];
        	while (ChildLeft(v) < (int)heap.size()){
            		int child = ChildRight(v) < (int)heap.size() && lt(heap[ChildRight(v)], heap[ChildLeft(v)]) ? ChildRight(v) : ChildLeft(v);
            		if (!lt(heap[child], x)) break;
            		heap[v] = heap[child];
			pos[heap[v]] = v;
			v = child;
        	}
        	heap[v] = x;
		pos[x] = v;
    	}

public:
    	void setComp   (Compare c)           { lt = c; }
    	bool empty     ()              const { return heap.size() == 0; }
    	bool inHeap    (int n)         const { return n < (int)pos.size() && pos[n] >= 0; }
    	void update    (int x)               { up(pos[x]); }

    	void insert(int x) {
        	if ((int)pos.size() < x + 1) pos.resize(x + 1, -1);
        	pos[x] = heap.size();
        	heap.push_back(x);
        	up(pos[x]); 
    	}

    	int pop() {
        	int x = heap[0];
        	heap[0] = heap.back();
        	pos[heap[0]] = 0;
		pos[x] = -1;
        	heap.pop_back();
        	if (heap.size() > 1) down(0);
        	return x; 
    	}
};

// Class for expressing clauses
class Clause {
public:
    	int lbd; // Literal Block Distance based on Glucose
    	std::vector<int> literals; // Literals in this clause
    	Clause( int sz ): lbd(0) { literals.resize(sz); }
    	int& operator [] ( int index ) { return literals[index]; }
};

// Class for expressing two-watched literals (lazy data structure)
class WL {
public:
    	int idx_clause; // Clause index in clause database
    	int blocker; // To fast guess whether a clause is already satisfied. 
    	WL(): idx_clause(0), blocker(0) {}
    	WL( int c, int b ): idx_clause(c), blocker(b) {}
};

// Function for soting activities
struct GreaterActivity { 
    	const double *activity;     
    	bool operator() ( int a, int b ) const { return activity[a] > activity[b]; }
    	GreaterActivity(): activity(NULL) {}
    	GreaterActivity( const double *s ): activity(s) {}
};


// Solver
class Solver {
public:
    	std::vector<int> learnt,                        // The clause indices of the learnt clauses.
                         trail,                         // Save the assigned literal sequence.
                         pos_in_trail,                  // Save the decision variables' position in trail.
                         reduce_map;                    // Auxiliary data structure for clause management.
    	std::vector<Clause> clause_DB;                  // clause database.
    	std::vector<WL> *watchedPointers;               // A mapping from literal to clauses.
    	int vars, clauses, origin_clauses, conflicts;   // the number of variables, clauses, conflicts.
    	int restarts, rephases, reduces;                // the number of conflicts since the last ... .
    	int rephase_limit, reduce_limit;                // parameters for when to conduct rephase and reduce.
    	int threshold;                                  // A threshold for updating the local_best phase.
    	int propagated;                                 // The number of propagted literals in trail.
    	int time_stamp;                                 // Aid parameter for conflict analyzation and LBD calculation.   
   
    	int lbd_queue[50],                              // circled queue saved the recent 50 LBDs.
            lbd_queue_size,                             // The number of LBDs in this queue
            lbd_queue_pos;                              // The position to save the next LBD.
    	double fast_lbd_sum, slow_lbd_sum;              // Sum of the Global and recent 50 LBDs.        
    	int *value,                                     // The variable assignement (1:True; -1:False; 0:Undefine) 
            *reason,                                    // The index of the clause that implies the variable assignment.
            *level,                                     // The decision level of a variable      
            *mark,                                      // Aid for conflict analyzation.
            *local_best,                                // A phase with a local deepest trail.                     
            *saved;                                     // Phase saving.
    	double *activity;                               // The variables' score for VSIDS.   
    	double var_inc;                                 // Parameter for VSIDS.               
    	Heap<GreaterActivity> vsids;                    // Heap to select variable.
     
    	void alloc_memory();                                      // Allocate memory for EasySAT 
    	void assign( int literal, int level, int cref );          // Assigned a variable.
    	int  propagate();                                         // BCP
    	void backtrack( int backtrack_level );                    // Backtracking
    	int  analyze( int cref, int &backtrack_level, int &lbd ); // Conflict analyzation.
    	int  parse( char *filename );                             // Read CNF file.
    	int  solve();                                             // Solving.
    	int  decide();                                            // Pick desicion variable.
    	int  add_clause( std::vector<int> &c );                   // add new clause to clause database.
    	void bump_var( int var, double mult );                    // update activity      
    	void restart();                                           // do restart.                                      
    	void reduce();                                            // do clause management.
    	void rephase();                                           // do rephase.
    	void printModel();                                        // print model when the result is SAT.
};
