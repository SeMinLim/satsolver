#include <stdio.h>
#include <sys/resource.h>

#include "solver.h"


// Elapsed time checker
static inline double timeCheckerCPU(void) {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}


// Main
int main() {
        Solver s;

        double parseTimeInitCPU = timeCheckerCPU();

        int res = solver_parse(&s);

        double parseTimeFinishCPU = timeCheckerCPU();

        double timeParseCPU = parseTimeFinishCPU - parseTimeInitCPU;

        printf( "Parsing Time (CPU): %.2f\n", timeParseCPU );

        if ( res == 20 ) printf("UNSATISFIABLE\n");
        else {
                double processTimeInitCPU = timeCheckerCPU();

                res = solver_solve(&s);
                if ( res == 10 ) {
			printf("SATISFIABLE\n");
                }
                else if ( res == 20 ) printf("UNSATISFIABLE\n");

                double processTimeFinishCPU = timeCheckerCPU();

                double timeProcessCPU = processTimeFinishCPU - processTimeInitCPU;

                printf( "Elapsed Time (CPU): %.2f\n", timeProcessCPU );
                printf( "Conflicts: %d\n", s.conflicts );
		printf( "Decisions: %d\n", s.decides );
		printf( "Propagations: %d\n", s.propagations );
		printf( "Evaluations: %d\n", s.propagations + s.decides );
        }
        return 0;
}
