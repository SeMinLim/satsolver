#ifndef HEAP_H
#define HEAP_H

#include <stdio.h>
#include <stdlib.h>

#include "vector.h"


#define ChildLeft(x) (x << 1 | 1)
#define ChildRight(x) ((x + 1) << 1)
#define Parent(x) ((x - 1) >> 1)


// Heap data structure
typedef struct GreaterActivity {
    	const double *activity;
    	int compare(int a, int b);
} GreaterActivity;

struct Heap {
    	Vector.init(heap);
    	Vector.init(pos);

    	void up( int v );
    	void down( int v );
	int empty();
	int inHeap(int n);
	void update(int x);
    	void insert(int x);
    	int pop();
} Heap;

#endif
