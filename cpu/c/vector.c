void VectorFunctions::init( Vector *v ) {
	v->capacity = VECTOR_INIT_CAPACITY;
	v->total = 0;
	v->items = malloc(sizeof(void*)*v->capacity);
}

int VectorFunctions::size( Vector *v ) {
	return v->total;
}

static void VectorFunctions::recapacity( Vector *v, int capacity ) {
	void **items = realloc(v->items, sizeof(void*)*capacity);
	if ( items ) {
		v->items = items;
		v->capacity = capacity;
	}
}

void VectorFunctions::resize( Vector *v, int size, int value ) {
	if ( size == v->total ) return;
	else if ( size > v->total ) while ( (size - v->total) != 0 ) pushback(v, value);
	else while ( (v->total - size) != 0 ) popback(v);
}
		
void VectorFunctions::pushback( Vector *v, void *item ) {
	if ( v->capacity == v->total ) recapacity(v, v->capacity*2);
	v->items[v->total++] = item;
}

void VectorFunctions::set( Vector *v, int index, void *item ) {
	if ( index >= 0 && index < v->total ) v->items[index] = item;
}

void VectorFunctions::*get( Vector *v, int index ) {
	if ( index >= 0 && index < v->total ) return v->items[index];
	return NULL;
}

void VectorFunctions::*back( Vector *v ) {
	return v->items[v->total];
}

void VectorFunctions::delete( Vector *v, int index ) {
	if ( index < 0 || index >= v->total ) return;
	
	v->items[index] = NULL;
	
	for ( int i = 0; i < v->total - 1; i++ ) {
		v->items[i] = v->items[i + 1];
		v->items[i + 1] = NULL;
	}
	
	v->total--;
	
	if ( v->total > 0 && v->total == v->capacity / 4 ) recapacity(v, v->capacity / 2);
}

void VectorFunctions::popback( Vector *v ) {
	v->items[v->total] = NULL;
	v->total--;

	if ( v->total > 0 && v->total == v->capacity / 4 ) recapacity(v, v->capacity / 2);
}

void VectorFunctions::clear( Vector *v ) {
	free(v->items);
}

