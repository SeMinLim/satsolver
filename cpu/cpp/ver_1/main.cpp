#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "solver.h"


// Elapsed time checker
static inline double timeCheckerCPU(void) {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

static inline double timeCheckerREAL() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return (double)tv.tv_sec + (double) tv.tv_usec / 1000000;
}


int main( int argc, char **argv ) {
        Solver S;

        double parseTimeInitCPU = timeCheckerCPU();
        double parseTimeInitREAL = timeCheckerREAL();

        int res = S.parse(argv[1]);

        double parseTimeFinishCPU = timeCheckerCPU();
        double parseTimeFinishREAL = timeCheckerREAL();

        double timeParseCPU = parseTimeFinishCPU - parseTimeInitCPU;
        double timeParseREAL = parseTimeFinishREAL - parseTimeInitREAL;

        printf( "Parsing Time (CPU): %.2f\n", timeParseCPU );
        printf( "Parsing Time (REAL): %.2f\n", timeParseREAL );

        if ( res == 20 ) printf("s UNSATISFIABLE\n");
        else {
                double processTimeInitCPU = timeCheckerCPU();
                double processTimeInitREAL = timeCheckerREAL();

                res = S.solve();
                if ( res == 10 ) {
			printf("s SATISFIABLE\n");
                        //S.printModel();
                }
                else if ( res == 20 ) printf("s UNSATISFIABLE\n");

                double processTimeFinishCPU = timeCheckerCPU();
                double processTimeFinishREAL = timeCheckerREAL();

                double timeProcessCPU = processTimeFinishCPU - processTimeInitCPU;
		double timeProcessREAL = processTimeFinishREAL - processTimeInitREAL;

                printf( "Elapsed Time (CPU): %.2f\n", timeProcessCPU );
                printf( "Elapsed Time (REAL): %.2f\n", timeProcessREAL );
		printf( "Conflicts: %d\n", S.conflicts );
		printf( "Decisions: %d\n", S.decides );
		printf( "Propagations: %d\n", S.propagations );
		printf( "Evaluations: %d\n", S.propagations + S.decides );
        }
        return 0;
}

