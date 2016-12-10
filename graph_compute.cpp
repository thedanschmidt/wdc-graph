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
int owner(long v, long n_proc, long num_vertices) {
    long own = (v*n_proc) / num_vertices;
    return (int) own;
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
    int bfs_idx = read_int( argc, argv, "-i", 0 );
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
    int max_vertex = (num_vertices / n_proc)*(rank+1)-1;
    if (max_vertex > num_vertices)
        max_vertex = num_vertices;
    
    printf("%d owns %d through %d \n", rank, min_vertex, max_vertex);

    // For now, iterate thorugh the whole edges file. If this is slow we can try a 
    // partitioned scheme
    std::ifstream infile(edge_file);
    if (infile == NULL) {
        printf("Failed to open edge file");
        exit(1);
    }

    int in_edge = 0;
    int out_edge = 0;
    int last_edge = min_vertex-1;
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
    row_ptr.push_back(col_index.size());
    
    // Check graph construction
    /*
    for (int i=0; i<col_index.size(); i++)
        printf("%d ", col_index[i]);
    printf("\n");
    for (int i=0; i<row_ptr.size(); i++)
        printf("%d ", row_ptr[i]);
    printf("\n");
    */

    double start_time = read_timer();

    if ( ! computation ) {
        printf("Computation option not specified.");
        exit(1);
    }
    else if (strcmp(computation, "bfs") == 0) {
        std::vector<int> distances(max_vertex-min_vertex+1, -1);
        int level = 1;
        std::stack<int> fs;
        std::stack<int> ns;
        if (n_proc == 1) {
            distances[bfs_idx] = 0;
            fs.push(bfs_idx);
            // Serial implementation of BFS
            while (!fs.empty()) {
                int u =0;
                while ( !fs.empty() ) {
                    u=fs.top();
                    for (int c = row_ptr[u]; c < row_ptr[u+1]; c++) {
                        int v = col_index[c];
                        if (distances[v] < 0) {
                            ns.push(v);
                            distances[v] = level;
                        }
                    }
                    fs.pop();
                }

                // TODO swap stacks with references for speed 
                fs = ns;
                while (!ns.empty()) {
                    ns.pop();
                }
                level++;
            }

            FILE* f = fopen("serial.out", "w");
            if (f != NULL) {
                for (int i=0; i<distances.size(); i++) {
                    fprintf(f, "%d\n", distances[i]);
                }
                fclose(f);
            }
            else {
                printf("failure to write serial out");
                exit(1);
            }
        }
        else {
            // Parallel Breadth First Search

            // Buffers for Alltoallv
            int* send_buf = new int[2*num_vertices];
            int* recv_buf = new int[2*num_vertices];
            int* send_counts = new int[n_proc];
            int* send_displs = new int[n_proc];
            int* recv_counts = new int[n_proc];
            int* recv_displs = new int[n_proc];

            // A data structure for making sure buffers are unique
            std::vector< std::set<int> > buff_sets;
            buff_sets.resize(n_proc);

            int* all_distances = NULL;
            if (rank == 0)
                all_distances = new int[num_vertices];

            for (int i=0; i<n_proc; i++) {
                send_displs[i] = i*(num_vertices/n_proc)+1;
                recv_displs[i] = i*(num_vertices/n_proc)+1;
                recv_counts[i] = 0;
            }
            int own = owner(bfs_idx, n_proc, num_vertices);
            if (own == rank) {
                fs.push(bfs_idx);
                distances[bfs_idx-min_vertex] = 0;
            }
            int not_finished = true;
            do {

                for (int i=0; i<n_proc; i++) {
                    send_counts[i] = 0;
                }
                while (!fs.empty()) {
                    int u=fs.top(); 
                    fs.pop();
                    for (int c = row_ptr[u-min_vertex]; c < row_ptr[u+1-min_vertex]; c++) {
                        int v = col_index[c];
                        int o = owner(v, n_proc, num_vertices);
                        if (o == rank) {
                            if (distances[v-min_vertex] < 0) {
                                distances[v-min_vertex] = level;
                                ns.push(v);
                            }
                        }
                        else {
                            // Put the vertex in the send buffer
                            if (buff_sets[o].find(v) != buff_sets[o].end()) {
                                int offset = send_displs[o]+send_counts[o];
                                send_buf[offset] = v;
                                send_counts[o]++;
                                buff_sets[o].insert(v);
                            }
                        }
                    }
                }

                // All-to-all sync fs
                // First exchange counts 
                MPI_Alltoall(
                    send_counts,
                    1,
                    MPI_INT,
                    recv_counts,
                    1,
                    MPI_INT,
                    MPI_COMM_WORLD
                );

                // Now exchange buffers
                MPI_Alltoallv(
                    send_buf,
                    send_counts,
                    send_displs,
                    MPI_INT,
                    recv_buf,
                    recv_counts,
                    recv_displs,
                    MPI_INT,
                    MPI_COMM_WORLD
                );

                // Search vertices from all-to-all 
                for (int p=0; p<n_proc; p++) {
                    if (p != rank) {
                        for (int r=0; r<recv_counts[p]; r++) {
                            int offset = recv_displs[p]+r;
                            int u = recv_buf[offset];
                            if (distances[u-min_vertex] < 0 ) {
                                distances[u-min_vertex] = level;
                                ns.push(u);
                            }
                        }
                    }
                }

                // TODO swap stacks with references for speed 
                while (!ns.empty()) {
                    fs.push(ns.top());
                    ns.pop();
                }
                level++;

                // Check if all proesses have finished
                int is_empty = fs.empty(); 
                MPI_Allreduce(
                    &(is_empty),
                    &not_finished,
                    1,
                    MPI_INT,
                    MPI_LAND,
                    MPI_COMM_WORLD
                );

            } while (!not_finished);
            printf("%d finished\n", rank);
            
            // Collect all distances on root process
            MPI_Gather(
                distances.data(),
                distances.size(),
                MPI_INT,
                all_distances,
                distances.size()+1,
                MPI_INT,
                0,
                MPI_COMM_WORLD
            );

            if (rank == 0) {
                FILE* f = fopen("parallel.out", "w");
                if (f != NULL) {
                    for (int i=0; i<num_vertices; i++) {
                        fprintf(f, "%d\n", all_distances[i]);
                    }
                    fclose(f);
                }
                else {
                    printf("failure to write parallel out");
                    exit(1);
                }
            }

            delete [] send_buf;
            delete [] recv_buf;
            delete [] send_counts;
            delete [] send_displs;
            delete [] recv_counts;
            delete [] recv_displs;
            delete [] all_distances;
        }
    }
    

    double end_time = read_timer();
    double time = end_time - start_time;
  
    if (rank == 0) {  
        // Summarize the results
        printf("Number of processes: %d took %f\n", n_proc, time);
        printf("Start time: %f, End time, %f\n", start_time, end_time);

    }
  
    //  release resources
    MPI_Finalize( );
    return 0;
}
