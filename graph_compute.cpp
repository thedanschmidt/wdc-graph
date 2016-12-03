#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <cstring>
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
        printf( "-nv Number of vertices to process on hyperlink graph \n" );
        printf( "-d Location of data\n" );
        return 0;
    }

    int num_vertices = read_int( argc, argv, "-nv", 1000 );
    const char* edge_file = read_string( argc, argv, "-d", NULL );

    std::vector<int> col_index;
    std::vector<int> row_ptr;
    
    //  set up MPI
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // Vertices owned by this process
    int min_vertex = (num_vertices / n_proc)*rank;
    int max_vertex = (num_vertices / n_proc)*(rank+1);
    if (max_vertex > num_vertices)
        max_vertex = num_vertices;
    
    printf("%d owns %d through %d \n", rank, min_vertex, max_vertex);

    // Global info everyone needs
    //MPI_Datatype PARTICLE;
    //MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    //MPI_Type_commit( &PARTICLE );
     

    // For now, iterate thorugh the whole edges file. If this is slow we can try a partitioned scheme
    FILE* ef = fopen(edge_file, "r");
    if (ef == NULL) {
        printf("Failed to open edge file");
        exit(1);
    }
    char* line;
    size_t len = 0;
    int out_edge = 0;
    int in_edge = 0;
    int last_edge = NULL; // is there a clver way to do this?
    int count = 0;

    while(getline(&line, &len, ef)) {
        in_edge = atoi(strtok(line, "\t")); // advances the line ptr
        if (in_edge >= min_vertex) {
            out_edge = atoi(line);
            if (in_edge > max_vertex) // bounds check
                break;

            printf("%d to %d\n", in_edge, out_edge);

            // Data structure construction
            // TODO: this doesn't work yet -- need to have a value for last_edge.
            if (in_edge == last_edge) {
                count++;
                col_index.push_back(out_edge);
            } else {
                row_ptr.push_back(row_ptr.back() + count);
                count = 0;
                col_index.push_back(out_edge);
            }

            // write 2 accessor functions --


        }
    }


    // TODO: serial BFS
    // input vertex index, return distance output indices out_vertex\tdistance
    // openMP shared queue...
    double time = read_timer( );

    // TODO Process graph here

    time = read_timer( ) - time;
  
    if (rank == 0) {  
        // Summarize the results
        printf("Number of processes: %d \n", n_proc);
    }
  
    //  release resources
    MPI_Finalize( );
    return 0;
}
