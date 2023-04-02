// greater_activity struct
bool greater_activity::compare( int a, int b ) {
	return activity[a] > activity[b];
}

// heap struct
void heap::up( int v ) {
        int x = VECTOR_GET(heap, int, v);
	int p = PARENT(v);
		
	// Child > Parent -> True
        while ( v && greater_activity.compare(x, VECTOR_GET(heap, int, p)) ) {
		VECTOR_SET(heap, v, VECTOR_GET(heap, int, p));
		VECTOR_SET(pos, VECTOR_GET(heap, int, p), v);
           	v = p;
		p = PARENT(p);
        }
        
	VECTOR_SET(heap, v, x);
	VECTOR_SET(pos, x, v);
}

void heap::down( int v ) {
	int x = VECTOR_GET(heap, int, v);

	while ( CHILDLEFT(v) < VECTOR_SIZE(heap) ) {
		// Pick the bigger one among left and right child
		int child = (CHILDRIGHT(v) < VECTOR_SIZE(heap)) && 
			    greater_activity.compare(VECTOR_GET(heap, int, CHILDRIGHT(v)), VECTOR_GET(heap, int, CHILDLEFT(v))) ? 
			    ChILDRIGHT(v) : CHILDLEFT(v);
            	if ( greater_activity.compare(x, VECTOR_GET(heap, int, child)) ) break;
		else {
			VECTOR_SET(heap, v, VECTOR_GET(heap, int, child));
			VECTOR_SET(pos, VECTOR_GET(heap, int, v), v);
			v = child;
		}
        }
	VECTOR_SET(heap, v, x);
	VECTOR_SET(pos, x, v);
}

bool heap::empty() { 
	return VECTOR_SIZE(heap) == 0; 
}

bool heap::inHeap( int n ) { 
	return (n < VECTOR_SIZE(pos)) && (VECTOR_GET(pos, int, n) >= 0); 
}

void heap::update( int x ) { 
	up(VECTOR_GET(pos, int, x)); 
}

void heap::insert( int x ) {
	if ( VECTOR_SIZE(pos) < x + 1 ) VECTOR_RESIZE(pos, x+1, -1);
	VECTOR_SET(pos, x, VECTOR_SIZE(heap));
	VECTOR_PUSHBACK(heap, x);
        up(VECTOR_GET(pos, int, x));
}

int heap::pop() {
	int x = VECTOR_GET(heap, int, 0);
	VECTOR_SET(heap, 0, VECTOR_BACK(heap));

	VECTOR_SET(pos, VECTOR_GET(heap, int, 0), 0);
	VECTOR_SET(pos, x, -1);
	VECTOR_POPBACK(heap);

        if ( VECTOR_SIZE(heap) > 1 ) down(0);
        return x;
}
