#ifndef VECTOR_H
#define VECTOR_H

#include <stdio.h>
#include <stdlib.h>


#define VECTOR_INIT_CAPACITY 4

#define Vector.init(vec) Vector vec; VectorFunctions.init(&vec)
#define Vector.size(vec) VectorFunctions.size(&vec)
#define Vector.resize(vec, size, value) VectorFunctions.resize(&vec, size, value)
#define Vector.pushback(vec, item) VectorFunctions.pushback(&vec, (void*) item)
#define Vector.set(vec, id, item) VectorFunctions.set(&vec, id, (void*) item)
#define Vector.get(vec, type, id) (type*) VectorFunctions.get(&vec, id)
#define Vector.back(vec, type) (type*) VectorFunctions.back(&vec)
#define Vector.delete(vec, id) VectorFunctions.delete(&vec, id)
#define Vector.popback(vec) VectorFunctions.popback(&vec)
#define Vector.clear(vec) VectorFunctions.clear(&vec)


typedef struct Vector {
	void **items;
	int capacity;
	int total;
} Vector;

typedef struct VectorFunctions {
	void init( Vector *v );
	int size( Vector *v );
	static void recapacity( Vector *v, int capacity );
	void resize( Vector *v, int size, int value );
	void pushback( Vector *v, void *item );
	void set( Vector *v, int index, void *item );
	void *get( Vector *v, int index );
	void *back( Vector *v );
	void delete( Vector *v, int index );
	void popback( Vector *v );
	void clear( Vector *v );
} VectorFunctions;

#endif
