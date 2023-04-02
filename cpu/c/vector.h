#ifndef VECTOR_H
#define VECTOR_H

#include <stdio.h>
#include <stdlib.h>


#define VECTOR_INIT_CAPACITY 4

#define VECTOR_INIT(vec) vector vec; vector_functions.init(&vec)
#define VECTOR_SIZE(vec) vector_functions.size(&vec)
#define VECTOR_RESIZE(vec, size, value) vector_functions.resize(&vec, size, value)
#define VECTOR_PUSHBACK(vec, item) vector_functions.pushback(&vec, (void*) item)
#define VECTOR_SET(vec, id, item) vector_functions.set(&vec, id, (void*) item)
#define VECTOR_GET(vec, type, id) (type*) vector_functions.get(&vec, id)
#define VECTOR_BACK(vec, type) (type*) vector_functions.back(&vec)
#define VECTOR_DELETE(vec, id) vector_functions.delete(&vec, id)
#define VECTOR_POPBACK(vec) vector_functions.popback(&vec)
#define VECTOR_CLEAR(vec) vector_functions.clear(&vec)


typedef struct vector {
	void **items;
	int capacity;
	int total;
} vector;

typedef struct vector_functions {
	void init( vector *v );
	int size( vector *v );
	static void recapacity( vector *v, int capacity );
	void resize( vector *v, int size, int value );
	void pushback( vector *v, void *item );
	void set( vector *v, int index, void *item );
	void *get( vector *v, int index );
	void *back( vector *v );
	void delete( vector *v, int index );
	void popback( vector *v );
	void clear( vector *v );
} vector_functions;

#endif
