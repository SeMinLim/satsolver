#include "solver.h"


char volatile * const printchar = (char*)0xfff;


int main( int argc, char **argv ) {
        Solver S;

        int res = S.parse(argv[1]);

        if ( res == 20 ) {
		//printf("s UNSATISFIABLE\n");
		(*printchar) = 0x55; // U
		(*printchar) = 0x4e; // N
		(*printchar) = 0x53; // S
		(*printchar) = 0x41; // A
		(*printchar) = 0x54; // T
		(*printchar) = 0x49; // I
		(*printchar) = 0x53; // S
		(*printchar) = 0x46; // F
		(*printchar) = 0x49; // I
		(*printchar) = 0x41; // A
		(*printchar) = 0x42; // B
		(*printchar) = 0x4c; // L
		(*printchar) = 0x45; // E
		(*printchar) = 0x0a; // \n
	}
        else {
                res = S.solve();
                if ( res == 10 ) {
                        //printf("s SATISFIABLE\n");
			(*printchar) = 0x53; // S
			(*printchar) = 0x41; // A
			(*printchar) = 0x54; // T
			(*printchar) = 0x49; // I
			(*printchar) = 0x53; // S
			(*printchar) = 0x46; // F
			(*printchar) = 0x49; // I
			(*printchar) = 0x41; // A
			(*printchar) = 0x42; // B
			(*printchar) = 0x4c; // L
			(*printchar) = 0x45; // E
			(*printchar) = 0x0a; // \n
                }
                else if ( res == 20 ) {
			//printf("s UNSATISFIABLE\n");
			(*printchar) = 0x55; // U
			(*printchar) = 0x4e; // N
			(*printchar) = 0x53; // S
			(*printchar) = 0x41; // A
			(*printchar) = 0x54; // T
			(*printchar) = 0x49; // I
			(*printchar) = 0x53; // S
			(*printchar) = 0x46; // F
			(*printchar) = 0x49; // I
			(*printchar) = 0x41; // A
			(*printchar) = 0x42; // B
			(*printchar) = 0x4c; // L
			(*printchar) = 0x45; // E
			(*printchar) = 0x0a; // \n
		}
        }
}

