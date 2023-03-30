#include <stdio.h>
#include <stdlib.h>


#define VECTOR_INIT_CAPACITY 4

#define VECTOR_INIT(vec) vector vec; vector_init(&vec)
#define VECTOR_SIZE(vec) vector_size(&vec)
#define VECTOR_PUSHBACK(vec, item) vector_pushback(&vec, (void *) item)
#define VECTOR_SET(vec, id, item) vector_set(&vec, id, (void *) item)
#define VECTOR_GET(vec, type, id) (type) vector_get(&vec, id)
#define VECTOR_DELETE(vec, id) vector_delete(&vec, id)
#define VECTOR_CLEAR(vec) vector_clear(&vec)


typedef struct vector {
	void **items;
	int capacity;
	int total;
} vector;

void vector_init( vector *v );
int vector_size( vector *v );
static void vector_resize( vector *v, int capacity );
void vector_pushback( vector *v, void *item );
void vector_set( vector *v, int index, void *item );
void *vector_get( vector *v, int index );
void vector_delete( vector *v, int index );
void vector_clear( vector *v );
