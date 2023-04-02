#ifndef HEAP_H
#define HEAP_H

#include <stdio.h>
#include <stdlib.h>

#include "vector.h"


#define CHILDLEFT(x) (x << 1 | 1)
#define CHILDRIGHT(x) ((x + 1) << 1)
#define PARENT(x) ((x - 1) >> 1)


// Heap data structure
typedef struct greater_activity {
    	const double *activity;
    	bool compare(int a, int b);
} greater_activity;

struct heap {
    	VECTOR_INIT(heap);
    	VECTOR_INIT(pos);

    	void up( int v );
    	void down( int v );
	bool empty();
	bool inHeap(int n);
	void update(int x);
    	void insert(int x);
    	int pop();
} heap;

#endif
