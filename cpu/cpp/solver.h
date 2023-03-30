#include <vector>
#include <fstream>

#define ChildLeft(x) (x << 1 | 1)
#define ChildRight(x) ((x + 1) << 1)
#define Parent(x) ((x - 1) >> 1)


// Heap data structure
template<class Compare>
class Heap {
    	Compare lt; // If left one of 'compare function' is big, return 'True'
    	std::vector<int> heap; // Value
    	std::vector<int> pos; // Position
    
    	void up(int v) {
        	int x = heap[v];
		int p = Parent(v);
		// Child > Parent -> True
        	while ( v && lt(x, heap[p]) ) {
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
            		// Pick the bigger one among left and right child
			int child = (ChildRight(v) < (int)heap.size()) && lt(heap[ChildRight(v)], heap[ChildLeft(v)]) ? ChildRight(v) : ChildLeft(v);
            		if ( lt(x, heap[child]) ) break;
			else {
				heap[v] = heap[child];
				pos[heap[v]] = v;
				v = child;
			}
        	}
        	heap[v] = x;
		pos[x] = v;
    	}

public:
    	void setComp   (Compare c)           { lt = c; }
    	bool empty     ()              const { return heap.size() == 0; }
    	bool inHeap    (int n)         const { return (n < (int)pos.size()) && (pos[n] >= 0); }
    	void update    (int x)               { up(pos[x]); }

    	void insert(int x) {
        	if ( (int)pos.size() < x + 1 ) pos.resize(x + 1, -1);	
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
        	if ( heap.size() > 1 ) down(0);
        	return x; 
    	}
};


// Clauses
class Clause {
public:
	// Literal block distance based on Glucose
	// LBD = How many decision variable in a learnt clause
    	int lbd;
    	std::vector<int> literals; // Literals in this clause
    	int& operator [] ( int index ) { return literals[index]; } // Overloading array operator
    	Clause( int sz ): lbd(0) { literals.resize(sz); } // Initialize int lbd as 0 and set a clause size with int size
};


// The watched literals data structure (lazy data structure)
class WL {
public:
    	int clauseIdx; // Clause index in clause database
    	int blocker; // To fast guess whether a clause is already satisfied. 
    	WL(): clauseIdx(0), blocker(0) {}
    	WL( int c, int b ): clauseIdx(c), blocker(b) {}
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
    	std::vector<int> learnt,                        // The clause indices of the learnt clauses
                         trail,                         // Save the assigned literal sequence(phase saving)
                         decVarInTrail,                 // Save the decision variables' position in trail(phase saving)
                         reduce_map;                    // Data structure for reduce
    	std::vector<Clause> clause_DB;                  // Clause database
    	std::vector<WL> *watchedPointers;               // A mapping from literal to clauses
    	
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
    	Heap<GreaterActivity> vsids;                    // Heap to select variable
     
    	void alloc_memory();                                      // Allocate memory 
    	void assign( int literal, int level, int cref );          // Assigned a variable
    	int  propagate();                                         // BCP (Boolean Contraint Propagation)
    	void backtrack( int backtrack_level );                    // Backtracking
    	int  analyze( int cref, int &backtrack_level, int &lbd ); // Conflict analyzation
    	int  parse( char *filename );                             // Read CNF file
    	int  solve();                                             // Solving
    	int  decide();                                            // Pick desicion variable
    	int  add_clause( std::vector<int> &c );                   // Add new clause to clause database
    	void bump_var( int var, double mult );                    // Update activity      
    	void restart();                                           // Do restart                                      
    	void reduce();                                            // Do reduce
    	void rephase();                                           // Do rephase
    	void printModel();                                        // Print model when the result is SAT
};
