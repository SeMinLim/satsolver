#include <time.h>
#include <sys/time.h>

#include "solver.h"


// Elapsed time checker
static inline double timeCheckerREAL() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return (double)tv.tv_sec + (double) tv.tv_usec / 1000000;
}


int main( int argc, char **argv ) {
        Solver S;

        double parseTimeInitREAL = timeCheckerREAL();

        int res = S.parse(argv[1]);

        double parseTimeFinishREAL = timeCheckerREAL();

        double timeParseREAL = parseTimeFinishREAL - parseTimeInitREAL;

        printf( "Parsing Time: %.2f\n", timeParseREAL );

        if ( res == 20 ) printf("s UNSATISFIABLE\n");
        else {
                double processTimeInitREAL = timeCheckerREAL();

                res = S.solve();
                if ( res == 10 ) {
                         printf("s SATISFIABLE\n");
                        S.printModel();
                }
                else if ( res == 20 ) printf("s UNSATISFIABLE\n");

                double processTimeFinishREAL = timeCheckerREAL();

                double timeProcessREAL = processTimeFinishREAL - processTimeInitREAL;

                printf( "Elapsed Time: %.2f\n", timeProcessREAL );
		printf( "Conflicts: %d\n", S.conflicts );
		printf( "Decisions: %d\n", S.decides );
		printf( "Propagations: %d\n", S.propagations );
		printf( "Evaluations: %d\n", S.propagations + S.decides );
        }
        return 0;
}

