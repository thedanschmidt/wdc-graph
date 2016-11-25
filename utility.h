#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include <vector>

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//  saving parameters
//
//  timing routines
//
double read_timer( );

//
//  I/O routines
//

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
