#include "solver.h"


//char volatile * const printchar = (char*)0x1fffffff;


// Main
int main() {
        Solver s;

        int res = solver_parse(&s);

        if ( res == 20 ) {
                //(*printchar) = 0x55; // U
                //(*printchar) = 0x4e; // N
                //(*printchar) = 0x53; // S
                //(*printchar) = 0x41; // A
                //(*printchar) = 0x54; // T
                //(*printchar) = 0xa; // \n
        } else {
                res = solver_solve(&s);
                if ( res == 10 ) {
                        //(*printchar) = 0x53; // S
                        //(*printchar) = 0x41; // A
                        //(*printchar) = 0x54; // T
                        //(*printchar) = 0xa; // \n
                } else if ( res == 20 ) {
                        //(*printchar) = 0x55; // U
                        //(*printchar) = 0x4e; // N
                        //(*printchar) = 0x53; // S
                        //(*printchar) = 0x41; // A
                        //(*printchar) = 0x54; // T
                        //(*printchar) = 0xa; // \n
                }
        }
}
