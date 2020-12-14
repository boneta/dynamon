#include <stdlib.h>
#include <time.h>
int my_random_ () { 
	srand( time( NULL ) ); 
	return( (int)( 32768. * rand() / ( RAND_MAX + 1. ) ) ); 
}
