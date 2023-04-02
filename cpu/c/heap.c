// Greater activity struct
int GreaterActivity::compare( int a, int b ) {
	if ( activity[a] > activity[b] ) return 1;
	else return 0;
}

// Heap struct
void Heap::up( int v ) {
        int x = Vector.get(heap, int, v);
	int p = Parent(v);
		
	// Child > Parent -> True
        while ( v && GreaterActivity.compare(x, Vector.get(heap, int, p)) ) {
		Vector.set(heap, v, Vector.get(heap, int, p));
		Vector.set(pos, Vector.get(heap, int, p), v);
           	v = p;
		p = Parent(p);
        }
        
	Vector.set(heap, v, x);
	Vector.set(pos, x, v);
}

void Heap::down( int v ) {
	int x = Vector.get(heap, int, v);

	while ( ChildLeft(v) < Vector.size(heap) ) {
		// Pick the bigger one among left and right child
		int child = (ChildRight(v) < Vector.size(heap)) && 
			    GreaterActivity.compare(Vector.get(heap, int, ChildRight(v)), Vector.get(heap, int, ChildLeft(v))) ? 
			    ChildRight(v) : ChildLeft(v);
            	if ( GreaterActivity.compare(x, Vector.get(heap, int, child)) ) break;
		else {
			Vector.set(heap, v, Vector.get(heap, int, child));
			Vector.set(pos, Vector.get(heap, int, v), v);
			v = child;
		}
        }
	Vector.set(heap, v, x);
	Vector.set(pos, x, v);
}

int Heap::empty() {
	if ( Vector.size(heap) == 0 ) return 1;
	else return 0;
}

int Heap::inHeap( int n ) { 
	if ( (n < Vector.size(pos)) && (Vector.get(pos, int, n) >= 0) ) return 1;
	else return 0;
}

void Heap::update( int x ) { 
	up(Vector.get(pos, int, x)); 
}

void Heap::insert( int x ) {
	if ( Vector.size(pos) < x + 1 ) Vector.resize(pos, x + 1, -1);
	Vector.set(pos, x, Vector.size(heap));
	Vector.pushback(heap, x);
        up(Vector.get(pos, int, x));
}

int Heap::pop() {
	int x = Vector.get(heap, int, 0);
	Vector.set(heap, 0, Vector.back(heap));

	Vector.set(pos, Vector.get(heap, int, 0), 0);
	Vector.set(pos, x, -1);
	Vector.popback(heap);

        if ( Vector.size(heap) > 1 ) down(0);
        return x;
}
