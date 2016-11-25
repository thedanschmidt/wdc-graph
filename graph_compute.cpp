#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <set>
#include "graph_compute.h"
#include "utility.h"

//
//  Compute basic properties of the hyperlink graph 
//
int main( int argc, char **argv )
{    
    //  process command line parameters
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        return 0;
    }
    
    //  set up MPI
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    // Global info everyone needs
    //MPI_Datatype PARTICLE;
    //MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    //MPI_Type_commit( &PARTICLE );

    double time = read_timer( );

    // TODO Process graph here

    time = read_timer( ) - time;
  
    if (rank == 0) {  
        // Summarize the results
        std::cout << "Number of processes: " << n_proc << std::endl;
    }
  
    //  release resources
    MPI_Finalize( );
    return 0;
}
