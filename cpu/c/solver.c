#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>


#define VECTOR_INIT_CAPACITY 4
#define VECTOR_INIT(vec) vector vec; vector_init(&vec)
#define VECTOR_PUSHBACK(vec, item) vector_pushback(&vec, (void *) item)
#define VECTOR_SET(vec, id, item) vector_set(&vec, id, (void *) item)
#define VECTOR_GET(vec, type, id) (type) vector_get(&vec, id)
#define VECTOR_DELETE(vec, id) vector_delete(&vec, id)
#define VECTOR_SIZE(vec) vector_size(&vec)
#define VECTOR_CLEAR(vec) vector_clear(&vec)


#define CHILDLEFT(x) (x << 1 | 1)
#define CHILDRIGHT(x) ((x + 1) << 1)
#define PARENT(x) ((x - 1) >> 1)


#define VALUE(literal) (literal > 0 ? value[literal] : -value[-literal])
#define WATCHPOINTER(id) (watchedPointers[vars + id])


// Functions for vector
typedef struct vector {
	void **items;
	int capacity;
	int total;
} vector;

void vector_init( vector *v ) {
	v->capacity = VECTOR_INIT_CAPACITY;
	v->total = 0;
	v->items = malloc(sizeof(void*)*v->capacity);
}

int vector_size( vector *v ) {
	return v->total;
}

static void vector_resize( vector *v, int capacity ) {
	void **items = realloc(v->items, sizeof(void*)*capacity);
	if ( items ) {
		v->items = items;
		v->capacity = capacity;
	}
}

void vector_pushback( vector *v, void *item ) {
	if ( v->capacity == v->total ) vector_resize(v, v->capacity * 2);
	v->items[v->total++] = item;
}

void vector_set( vector *v, int index, void *item ) {
	if ( index >= 0 && index < v->total ) v->items[index] = item;
}

void *vector_get( vector *v, int index ) {
	if ( index >= 0 && index < v->total ) return v->items[index];
	return NULL;
}

void vector_delete( vector *v, int index ) {
	if ( index < 0 || index >= v->total ) return;
	
	v->items[index] = NULL;
	
	for ( int i = 0; i < v->total - 1; i++ ) {
		v->items[i] = v->items[i + 1];
		v->items[i + 1] = NULL;
	}
	
	v->total--;
	
	if ( v->total > 0 && v->total == v->capacity / 4 ) vector_resize(v, v->capacity / 2);
}

void vector_clear( vector *v ) {
	free(v->items);
}


// Heap data structure
typedef struct greater_activity { 
    	const double *activity;     
    	bool operator() ( int a, int b ) const { return activity[a] > activity[b]; }
    	greaterActivity(): activity(NULL) {}
    	greaterActivity( const double *s ): activity(s) {}
} greater_activity;

typedef struct heap {
    	Compare lt; // TODO: fix this! // If left one of 'compare function' is big, return 'True'
    	VECTOR_INIT(heap); // Value
    	VECTOR_INIT(pos); // Position

    	void up( int v ) {
        	int x = VECTOR_GET(heap, v);
		int p = PARENT(v);
		
		// Child > Parent -> True
        	while ( v && lt(x, VECTOR_GET(heap, p)) ) {
			VECTOR_SET(heap, v, VECTOR_GET(heap, p));
			VECTOR_SET(pos, VECTOR_GET(heap, p), v);
            		v = p;
			p = PARENT(p);
        	}
        	
		VECTOR_SET(heap, v, x);
		VECTOR_SET(pos, x, v);
    	}

    	void down( int v ) {
		int x = VECTOR_GET(heap, v);
        	
		while ( CHILDLEFT(v) < VECTOR_SIZE(heap) ) {
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
} heap;


// Functions for parser
char *read_whitespace( char *p ) {
        // ASCII
        // Horizontal tab, line feed or new line, vertical tab, form feed or new page, carriage return, space
        while ( (*p >= 9 && *p <= 13) || *p == 32 ) ++p;
        return p;
}
char *read_until_new_line( char *p ) {
        while ( *p != '\n' ) {
                if ( *p++ == '\0' ) exit(1);
        }
        return ++p;
}
char *read_int( char *p, int *i ) {
        bool sym = true;
        *i = 0;
        p = read_whitespace(p);
        if ( *p == '-' ) {
                sym = false;
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


// Functions for Solver
int add_clause( vector &c ) {                   
    	clause_DB.push_back(Clause(c.size()));                          
    	int id = clause_DB.size() - 1;                                
    	for ( int i = 0; i < (int)c.size(); i++ ) clause_DB[id][i] = c[i];
        // There's two watched literals	
    	WatchPointer(-c[0]).push_back(WL(id, c[1])); // watchedPointers[vars-c[0]]                      
    	WatchPointer(-c[1]).push_back(WL(id, c[0])); // watchedPointers[vars-c[1]]
    	return id;                                                      
}

int parse( char *filename ) {
	FILE *f_data = fopen(filename, "r");
	
	fseek(f_data, 0, SEEK_END);
	size_t file_len = ftell(f_data); // Get the file size

	fseek(f_data, 0, SEEK_SET);
	char *data = (char*)malloc(sizeof(char)*(file_len+1));
	char *p = data;
	fread(data, sizeof(char), file_len, f_data);
	fclose(f_data);
	data[file_len] = '\0'; // Read the file
 
	VECTOR_INIT(buffer); // Save the clauses temporarily
	
    	while ( *p != '\0' ) {
        	p = read_whitespace(p);
        	if ( *p == '\0' ) break;
        	if ( *p == 'c' ) p = read_until_new_line(p); // If there are some comments in CNF
       	 	else if ( *p == 'p' ) { 
            		if ( (*(p + 1) == ' ') && (*(p + 2) == 'c') && 
			     (*(p + 3) == 'n') && (*(p + 4) == 'f') ) {
                		p += 5; 
				//p = read_int(p, &vars); 
				//p = read_int(p, &clauses);
                		//alloc_memory();
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
				else VECTOR_PUSHBACK(buffer, dimacs_lit);
			}
			else {                                                       
                		if ( VECTOR_SIZE(buffer) == 0 ) return 20;
				/*else if ( VECTOR_SIZE(buffer) == 1 ) {
					if ( Value(buffer[0]) == -1 ) return 20;
					else if ( !Value(buffer[0]) ) assign(buffer[0], 0, -1);
				}
                		else add_clause(buffer); 
				*/
                		//VECTOR_CLEAR(buffer);
   				//free(data);				
            		}
        	}
    	}
    	//origin_clauses = clause_DB.size();
    	return 10;//( propagate() == -1 ? 0 : 20 );             
}


// Elapsed time checker
static inline double timeCheckerCPU(void) {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

static inline double timeCheckerREAL() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return (double)tv.tv_sec + (double) tv.tv_usec / 1000000;
}

int main( int argc, char **argv ) {
        double parseTimeInitCPU = timeCheckerCPU();
        double parseTimeInitREAL = timeCheckerREAL();

        int res = parse(argv[1]);

        double parseTimeFinishCPU = timeCheckerCPU();
        double parseTimeFinishREAL = timeCheckerREAL();

        double timeParseCPU = parseTimeFinishCPU - parseTimeInitCPU;
        double timeParseREAL = parseTimeFinishREAL - parseTimeInitREAL;

        printf( "Parsing Time (CPU): %.2f\n", timeParseCPU );
        printf( "Parsing Time (REAL): %.2f\n", timeParseREAL );
	
/*	if ( res == 20 ) printf("s UNSATISFIABLE\n");
        else {
                double processTimeInitCPU = timeCheckerCPU();
                double processTimeInitREAL = timeCheckerREAL();

                res = S.solve();
                if ( res == 10 ) {
                         printf("s SATISFIABLE\n");
                        S.printModel();
                }
                else if ( res == 20 ) printf("s UNSATISFIABLE\n");

                double processTimeFinishCPU = timeCheckerCPU();
                double processTimeFinishREAL = timeCheckerREAL();

                double timeProcessCPU = processTimeFinishCPU - processTimeInitCPU;
                double timeProcessREAL = processTimeFinishREAL - processTimeInitREAL;

                printf( "Elapsed Time (CPU): %.2f\n", timeProcessCPU );
                printf( "Elapsed Time (REAL): %.2f\n", timeProcessREAL );
		printf( "Conflicts: %d\n", S.conflicts );
		printf( "Decisions: %d\n", S.decides );
		printf( "Propagations: %d\n", S.propagations );
		printf( "Evaluations: %d\n", S.propagations + S.decides );
        }
	*/
        return 0;
}
