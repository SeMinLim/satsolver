void vector_functions::init( vector *v ) {
	v->capacity = VECTOR_INIT_CAPACITY;
	v->total = 0;
	v->items = malloc(sizeof(void*)*v->capacity);
}

int vector_functions::size( vector *v ) {
	return v->total;
}

static void vector_functions::recapacity( vector *v, int capacity ) {
	void **items = realloc(v->items, sizeof(void*)*capacity);
	if ( items ) {
		v->items = items;
		v->capacity = capacity;
	}
}

void vector_functions::resize( vector *v, int size, int value ) {
	if ( size == v->total ) return;
	else if ( size > v->total ) while ( (size - v->total) != 0 ) pushback(v, value);
	else while ( (v->total - size) != 0 ) popback(v);
}
		
void vector_functions::pushback( vector *v, void *item ) {
	if ( v->capacity == v->total ) recapacity(v, v->capacity*2);
	v->items[v->total++] = item;
}

void vector_functions::set( vector *v, int index, void *item ) {
	if ( index >= 0 && index < v->total ) v->items[index] = item;
}

void vector_functions::*get( vector *v, int index ) {
	if ( index >= 0 && index < v->total ) return v->items[index];
	return NULL;
}

void vector_functions::*back( vector *v ) {
	return v->items[v->total];
}

void vector_functions::delete( vector *v, int index ) {
	if ( index < 0 || index >= v->total ) return;
	
	v->items[index] = NULL;
	
	for ( int i = 0; i < v->total - 1; i++ ) {
		v->items[i] = v->items[i + 1];
		v->items[i + 1] = NULL;
	}
	
	v->total--;
	
	if ( v->total > 0 && v->total == v->capacity / 4 ) recapacity(v, v->capacity / 2);
}

void vector_functions::popback( vector *v ) {
	v->items[v->total] = NULL;
	v->total--;

	if ( v->total > 0 && v->total == v->capacity / 4 ) recapacity(v, v->capacity / 2);
}

void vector_functions::clear( vector *v ) {
	free(v->items);
}

