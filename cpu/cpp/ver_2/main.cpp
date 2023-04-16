#include "solver.h"


int main( int argc, char **argv ) {
        Solver S;
        
        int res = S.parse(argv[1]);
        
        if ( res == 20 ) printf("UNSATISFIABLE\n");
        else {
                res = S.solve();
                if ( res == 10 ) {
			printf("SATISFIABLE\n");
                        //S.printModel();
                }
                else if ( res == 20 ) printf("UNSATISFIABLE\n");
		else if ( res == 30 ) printf( "UNSOLVED\n" );
		//printf( "Conflicts: %d\n", S.conflicts );
		//printf( "Decisions: %d\n", S.decides );
		//printf( "Propagations: %d\n", S.propagations );
		//printf( "Evaluations: %d\n", S.propagations + S.decides );
        }
        return 0;
}
