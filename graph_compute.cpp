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
#include <fstream>


/*
*
* Simple Accessor Function used in BFS
* Input: a vertex
* Output: a vector of vertices that are connected to the input vertex
*
*/
std::vector<int> findAllConnections(int in_vertex, std::vector<int>& row_ptr, std::vector<int>& col_index) {
    std::vector<int> connections;
    int startingIndex = row_ptr[in_vertex];
    int endingIndex = connections.size();
    if (in_vertex < (connections.size() - 1)) 
        endingIndex = row_ptr[in_vertex + 1];
    for (int i = startingIndex; i < endingIndex; i++) {
        connections.push_back(col_index[i]);
    }

    return connections;

}

// Gives the process that owns v
int owner(int v, int n_proc, int num_vertices) {
    int v_per = num_vertices / n_proc;
    return v / v_per;
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
     

    // For now, iterate thorugh the whole edges file. If this is slow we can try a 
    // partitioned scheme
    //FILE* ef = fopen(edge_file, "r");
    std::ifstream infile(edge_file);
    if (infile == NULL) {
        printf("Failed to open edge file");
        exit(1);
    }

    int in_edge = 0;
    int out_edge = 0;
    int last_edge = -1;
    int count = 0;

    while(infile >> in_edge >> out_edge) {
        //in_edge = atoi(strtok(line, "\t")); // advances the line ptr
        if (in_edge >= min_vertex) {
            if (in_edge > max_vertex)
                break;

            // Data structure construction
            if (out_edge <= num_vertices) {
                col_index.push_back(out_edge);
                if (in_edge != last_edge) {
                    while (last_edge != in_edge) {
                        row_ptr.push_back(count);
                        last_edge++;
                    }
                }         
                count++;
            }
        }
    }
    
    // Check graph construction
    /*
    for (int i=0; i<col_index.size(); i++)
        printf("%d ", col_index[i]);
    printf("\n");
    for (int i=0; i<row_ptr.size(); i++)
        printf("%d ", row_ptr[i]);
    printf("\n");
    */

    double time = read_timer( );

    if ( ! computation ) {
        printf("Computation option not specified.");
        exit(1);
    }
    else if (strcmp(computation, "bfs") == 0) {
        std::vector<int> distances(max_vertex-min_vertex, -1);
        int level = 1;
        std::stack<int> fs;
        std::stack<int> ns;
        if (bfs_idx >= min_vertex && bfs_idx <= max_vertex)
            distances[bfs_idx] = 0;
            fs.push(bfs_idx);
        if (n_proc == 1) {
            // Serial implementation of BFS
            while (!fs.empty()) {
                for (int u=fs.top(); !fs.empty(); fs.pop()) {
                    for (int c = row_ptr[u]; c < row_ptr[u+1]; c++) {
                        int v = col_index[c];
                        if (distances[v] < 0) {
                            ns.push(v);
                            distances[v] = level;
                        }
                    }
                }

                // TODO swap stacks with references for speed 
                while (!ns.empty()) {
                    fs.push(ns.top());
                    ns.pop();
                }
                level++;
            }
            for (int i=0; i<distances.size(); i++) {
                printf("%d\n", distances[i]);
            }
        }
        else {
            // Parallel Breadth First Search
            int* send_buf = new int[num_vertices];
            int* recv_buf = new int[num_vertices];
            do {
                for (int u=fs.top(); !fs.empty(); fs.pop()) {
                    for (int c = row_ptr[u]; c < row_ptr[u+1]; c++) {
                        int v = col_index[c];
                        int o = owner(v, n_proc, num_vertices);
                        if (o == rank) {
                            if (distances[v] < 0) {
                                ns.push(v);
                                distances[v] = level;
                            }
                        }
                        else {
                            // push on buffer
                        }
                    }
                    // All-to-all sync fs

                    // Search fs adding to ns
                   
                }
            } while (!fs.empty());

            delete [] send_buf;
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
