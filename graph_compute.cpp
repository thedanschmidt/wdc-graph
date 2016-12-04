#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <set>
#include <stack>
#include "graph_compute.h"
#include "utility.h"


/*
*
* Simple Accessor Function used in BFS
* Input: a vertex
* Output: a vector of vertices that are connected to the input vertex
*
*/
std::vector<int> findAllConnections(int in_vertex, std::vector<int> row_ptr, std::vector<int> col_index) {
    std::vector<int> connections = {};
    int startingIndex = row_ptr[in_vertex];
    int endingIndex = coonections.size();
    if (in_vertex < (connections.size() - 1)) 
        endingIndex = row_ptr[in_vertex + 1];
    for (int i = startingIndex; i < endingIndex; i++) {
        connections.push_back(col_ind[i]);
    }

    return connections;

}

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
        printf( "-c bfs|scc the computation to be performed");
        printf( "-i Start index for BFS");
        printf( "-o Output file");
        return 0;
    }

    int num_vertices = read_int( argc, argv, "-nv", 1000 );
    const char* computation = read_string(argc, argv, "-c", NULL);
    int bfs_idx = read_int( argc, argv, "-i", 1000 );
    const char* edge_file = read_string( argc, argv, "-d", NULL );

    std::vector<int> col_index;
    std::vector<int> row_ptr = {0};
    
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
     

    // For now, iterate thorugh the whole edges file. If this is slow we can try a 
    // partitioned scheme
    FILE* ef = fopen(edge_file, "r");
    if (ef == NULL) {
        printf("Failed to open edge file");
        exit(1);
    }
    char* line;
    size_t len = 0;
    int out_edge = 0;
    int in_edge = 0;
    int last_edge = -1;
    int count = 0;

    while(getline(&line, &len, ef)) {
        in_edge = atoi(strtok(line, "\t")); // advances the line ptr
        if (in_edge >= min_vertex) {
            if (last_edge == -1) 
                last_edge = in_edge;
            out_edge = atoi(line);
            if (in_edge > max_vertex)
                break;

            printf("%d to %d\n", in_edge, out_edge);

            // Data structure construction
            col_index.push_back(out_edge);
            if (in_edge == last_edge) {
                count++;
            } else {
                row_ptr.push_back(row_ptr.back() + count);
                count = 0;
            }

        if (count > 0) {
            row_ptr.push_back(row_ptr.back() + count);
            count = 0;
        }


        }
    }


    // TODO: serial BFS
    // input vertex index, return distance output indices out_vertex\tdistance
    // openMP shared queue...
    double time = read_timer( );

    std::vector<int> col_ind; 
    std::vector<int> row_ptr; 
    if ( ! computation ) {
        printf("Computation option not specified.");
        exit(1);
    }
    else if (strcmp(computation, "bfs") == 0) {
        if (n_proc == 1) {
            // Rohin code here :D
        }
        else {
            // Parallel Breadth First Search
            std::vector<int> distances;  
            distances.assign(max_vertex-min_vertex, num_vertices+1);
            distances[bfs_idx] = 0;
            int level = 1;
            std::stack<int> fs;
            std::stack<int> ns;
        }
    }
    

    time = read_timer( ) - time;
  
    if (rank == 0) {  
        // Summarize the results
        printf("Number of processes: %d \n", n_proc);
    }
  
    //  release resources
    MPI_Finalize( );
    return 0;
}
