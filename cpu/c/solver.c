#include <stdio.h>
#include <time.h>

#define MAXNUMLIT 5

double timespec_diff_sec( timespec start, timespec end ) {
	double t = end.tv_sec - start.tv_sec;
	t += ((double)(end.tv_nsec - start.tv_nsec)/1000000000L);
	return t;
}

// Functions for parser
char *read_whitespace( char *p ) {
        // ASCII
        // Horizontal tab, line feed or new line, vertical tab, form feed or new page, carriage return, space
        while ( (*p >= 9 && *p <= 13) || *p == 32 ) ++p;
        return p;
}
char *read_until_new_line( char *p ) {
        while ( *p != '\n' ) {
                if ( *p++ == '\0' ) exit(1);
        }
        return ++p;
}
char *read_int( char *p, int *i ) {
        bool sym = true;
        *i = 0;
        p = read_whitespace(p);
        if ( *p == '-' ) {
                sym = false;
                ++p;
        }
        while ( *p >= '0' && *p <= '9' ) {
                if ( *p == '\0' ) return p;
                *i = *i * 10 + *p - '0';
                ++p;
        }
        if ( !sym ) *i = -(*i);
        return p;
}


// Functions for Solver
void parse( char *filename  ) {
	FILE *f_data = fopen(filename, "r");
	
	fseek(f_data, 0, SEEK_END);
	size_t fileLength = ftell(f_data); // Get the file size

	fseek(f_data, 0, SEEK_SET);
	char *data = (char*)malloc(sizeof(char)*(fileLength+1));
	char *p = data;
	fread(data, sizeof(char), fileLength, f_data);
	fclose(f_data);
	data[fileLength] = '\0'; // Read the file

	// Save the clauses temporarily
	while( *p != '\0' ) {
	       	p = read_whitespace(p);
        	if ( *p == '\0' ) break;
        	if ( *p == 'c' ) p = read_until_new_line(p); // If there are some comments in CNF
       	 	else if ( *p == 'p' ) {
            		if ( (*(p + 1) == ' ') && (*(p + 2) == 'c') &&
			     (*(p + 3) == 'n') && (*(p + 4) == 'f') ) {
                		p += 5;
				p = read_int(p, &vars);
				printf("%d\n", p);
				p = read_int(p, &clauses);
                		//alloc_memory();
            		}
            		else printf("PARSE ERROR(Unexpected Char)!\n"), exit(2);
        	}
        	else {
            		int32_t dimacs_lit;
            		p = read_int(p, &dimacs_lit);
            		if ( dimacs_lit != 0 ) {
				if ( *p == '\0' ) {
                			printf("c PARSE ERROR(Unexpected EOF)!\n");
					exit(1);
				}
				//else buffer.push_back(dimacs_lit);
			}
			else {
                		/*if ( buffer.size() == 0 ) return 20;
				else if ( buffer.size() == 1 ) {
					if ( Value(buffer[0]) == -1 ) return 20;
					else if ( !Value(buffer[0]) ) assign(buffer[0], 0, -1);
				}
                		else add_clause(buffer);
				*/
                		//buffer.clear();
				return 20;
            		}
        	}
    	}
    	//origin_clauses = clause_DB.size();
    	//return ( propagate() == -1 ? 0 : 20 );
}

int main
