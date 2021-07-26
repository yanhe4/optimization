#include <iostream>
#include <fstream> // std::ifstream
#include <sstream> // std::stringstream
#include <string> // std::string, std::stoi
#include <cstring> // std::strcmp
#include <cmath>
#include <vector>
#include <chrono>
#include <ctime>
#include "Graph.hpp"
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/count.h>

#define GPU 1

using namespace optPagerank::gpu_graph;

extern const double damping_factor = 0.85;
extern const unsigned max_iterations = 100;
extern const double tolerance = 1e-10;
const int blocksize = 512;

// Read Input (pairs of source and destination links) from file with format:
// src_index dest_index
// ... 
// src_index dest_index 
ColdEdge ReadInputFromTextFile(const char* input_file, unsigned& num_vertices)
{
    std::ifstream myfile (input_file);
    ColdEdge edges;
    unsigned source, destination;
    if (myfile.is_open()) 
    {
      while(myfile >> source >> destination)
      {
        unsigned larger = (source > destination)? source : destination;
        num_vertices = (num_vertices > larger)? num_vertices : larger;
        edges.src.push_back(source);
        edges.dest.push_back(destination);  
      }
      ++num_vertices;
      myfile.close();
    }
    return edges;
}

bool ToleranceCheck(const unsigned& num_v, HotData& hotData)
{
    // Sum up the pagerank
    double pr_sum = 0.0;
    for (unsigned i = 0; i < num_v; i++) 
    {
        pr_sum += hotData.pagerank[i];
    }
    // Calculate the cur_toleranceor
    pr_sum = 1.0 / pr_sum;
    double cur_tolerance = 0.0;
    for (unsigned i = 0; i < num_v; i++)
    {
        hotData.pagerank[i] *= pr_sum;
        // norm 1
        cur_tolerance += std::fabs(hotData.pagerank[i] - hotData.pre_pagerank[i]);
    }
    if (cur_tolerance < tolerance)
    {
        std::cout << "Current toleranceor: " << cur_tolerance << std::endl;
        return true;
    }
    return false;
}

#ifdef GPU
__global__
void update_pagerank( int *ingoing_edges_num, int *outgoing_edges_num,
                      int *begin_index, int *adj_edges, double *pre_pagerank,
                      double pr_dangling, double* pr_random, double *pagerank,
                      double* dangling_vec, size_t n )
{
   int const index = threadIdx.x + blockIdx.x * blockDim.x;
   if( index < n ) {
      int num_edges = ingoing_edges_num[index];
      int begin_index_ = begin_index[index];
      pagerank[index] = 0.0;
      for( int i = 0; i < num_edges; ++i ){
         int inward_edge_index = adj_edges[begin_index_ + i];
         double pr_eigenvector = 0.85 * pre_pagerank[inward_edge_index]
                                 / outgoing_edges_num[inward_edge_index];
         pagerank[index] += pr_eigenvector;
      }
      pagerank[index] += (*pr_random + pr_dangling);
      dangling_vec[index] = pagerank[index] * (outgoing_edges_num[index] == 0);
   }
}
#endif

void PageRank(GPU_Graph *graph)
{
   const unsigned num_v = graph->VertexesNum();
   double init_rank = double(1.0 / num_v);
   double pr_random = (1.0 - damping_factor) / num_v;
   unsigned iter = 0;

#ifdef GPU

   // calculate number of blocks. block_size is fixed to 512
   auto const num_blocks = std::ceil( num_v / static_cast< float >( blocksize) );
   //Initialize all memories used by GPU
   int *dev_ingoing_edge_nums;
   int *dev_outgoing_edge_nums;
   int *dev_begin_index;
   int *dev_adj_edges;
   double *dev_pre_pagerank;
   double *dev_pagerank;
   double *dev_pr_random;
   double *dev_dangling_vec;

   //Allocate memory for elments mapped to GPU
   cudaMalloc( (void **) &dev_ingoing_edge_nums, num_v*sizeof(int) );
   cudaMalloc( (void **) &dev_outgoing_edge_nums, num_v*sizeof(int) );
   cudaMalloc( (void **) &dev_begin_index, num_v*sizeof(int) );
   cudaMalloc( (void **) &dev_adj_edges, graph->num_edges*sizeof(int) );
   cudaMalloc( (void **) &dev_pagerank, num_v*sizeof(double) );
   cudaMalloc( (void **) &dev_pre_pagerank, num_v*sizeof(double) );
   cudaMalloc( (void **) &dev_pr_random, sizeof(double) );
   cudaMalloc( (void **) &dev_dangling_vec, num_v*sizeof(double) );

   //Define trust objects and wrap raw pointers with device pointers
   thrust::device_ptr<double> dev_thrust_dangling_vec(dev_dangling_vec);

   //Initialize obejcts that won't be changed in the algorithm
   cudaMemcpy( dev_ingoing_edge_nums, graph->ingoing_edges_num.data(),
               num_v*sizeof(int), cudaMemcpyHostToDevice );
   cudaMemcpy( dev_outgoing_edge_nums, graph->hotData.outgoing_edges_num.data(),
               num_v*sizeof(int), cudaMemcpyHostToDevice );
   cudaMemcpy( dev_begin_index, graph->beginIndex.data(),
               num_v*sizeof(int), cudaMemcpyHostToDevice );
   cudaMemcpy( dev_adj_edges, graph->adjE,
               graph->num_edges*sizeof(int), cudaMemcpyHostToDevice );
   cudaMemcpy( dev_pr_random, &pr_random, sizeof(double), cudaMemcpyHostToDevice );

   //Initialize objects that will be updated in the algorithm, and copy
   //them from host to device
   graph->hotData.pagerank.assign(num_v, init_rank);
   graph->hotData.pre_pagerank.assign(num_v, 0);
   cudaMemcpy( dev_pre_pagerank, graph->hotData.pre_pagerank.data(),
               num_v*sizeof(double), cudaMemcpyHostToDevice );
   cudaMemcpy( dev_pagerank, graph->hotData.pagerank.data(),
               num_v*sizeof(double), cudaMemcpyHostToDevice );

   double dangling_pr_sum = 0.0;
   for (unsigned i = 0; i < num_v; i++)
   {
      dangling_pr_sum += graph->hotData.pagerank[i] * (graph->hotData.outgoing_edges_num[i] == 0);
   }

   while(iter++ < max_iterations){
      double pr_dangling = damping_factor * dangling_pr_sum / num_v;
      cudaMemcpy( dev_pre_pagerank, dev_pagerank, num_v*sizeof(double), cudaMemcpyDeviceToDevice );
      //Main function in this algorithm to update pagerank at each iteration, hand over to GPU
      update_pagerank<<< num_blocks, blocksize >>>(dev_ingoing_edge_nums, dev_outgoing_edge_nums,
                                                   dev_begin_index, dev_adj_edges, dev_pre_pagerank,
                                                   pr_dangling, dev_pr_random, dev_pagerank,
                                                   dev_dangling_vec, num_v);
      dangling_pr_sum = thrust::reduce( dev_thrust_dangling_vec, dev_thrust_dangling_vec+num_v );
   }
   cudaMemcpy( graph->hotData.pagerank.data(), dev_pagerank, num_v*sizeof(double), cudaMemcpyDeviceToHost );
   cudaFree( dev_ingoing_edge_nums );
   cudaFree( dev_outgoing_edge_nums );
   cudaFree( dev_begin_index );
   cudaFree( dev_adj_edges );
   cudaFree( dev_pagerank );
   cudaFree( dev_pre_pagerank );
   cudaFree( dev_pr_random );
   cudaFree( dev_dangling_vec );
#endif

#if 0
    // This is the original algorithm in CPU which we port to GPU
    for (unsigned i = 0; i < num_v; i++)
    {
        graph->hotData.pagerank[i] = init_rank;
        graph->hotData.pre_pagerank[i] = 0.0;
    }

    while (iter++ < max_iterations)
    {
        double dangling_pr_sum = 0.0;
        // Update the pagerank values in every iteration
        for (unsigned i = 0; i < num_v; i++)
        {
            graph->hotData.pre_pagerank[i] = graph->hotData.pagerank[i];
            graph->hotData.pagerank[i] = 0.0;
            dangling_pr_sum += graph->hotData.pre_pagerank[i] * (graph->hotData.outgoing_edges_num[i] == 0);
        }

        double pr_dangling = damping_factor * dangling_pr_sum / num_v;

        // Iterater all the vertexes and calculate its adjacency function l(pi,pj) of all inward links
        // Update its pagerank value by adding pr_eigenvector from its inward links separately
        for( int i = 0; i < num_v; ++i )
        {
            unsigned inward_edges_num = graph->ingoing_edges_num[i];
            int begin_index = graph->beginIndex[i];
            for( int j = 0; j < inward_edges_num; ++j){
               unsigned inward_edge_index = graph->adjE[begin_index + j];
               double pr_eigenvector = damping_factor * graph->hotData.pre_pagerank[inward_edge_index]
                                        / graph->hotData.outgoing_edges_num[inward_edge_index];
               graph->hotData.pagerank[i] += pr_eigenvector;
            }
            graph->hotData.pagerank[i] += (pr_random + pr_dangling);
        }
        // finish when cur_toleranceor is smaller than tolerance we set
        if(ToleranceCheck(num_v, graph->hotData)) 
        {
            std::cout << "Iteration time: " << iter << std::endl;
            break;
        }
    }
#endif
}

void printFinalResults(GPU_Graph* graph)
{
    std::cout << "PageRank values: \n";
    for(int i = 0; i < graph->VertexesNum(); ++i)
    {
        std::cout << "The index is: " << i << " with value " << graph->hotData.pagerank[i] << '\n';  
    }
    std::cout<<'\n';
}

void PrintBenchmark(std::chrono::time_point<std::chrono::steady_clock> start_t, std::chrono::time_point<std::chrono::steady_clock> const end_t, const unsigned loop_t)
{
    auto const avg_time = std::chrono::duration_cast<std::chrono::microseconds>( end_t - start_t ).count() / double(loop_t);
    std::cout << "Average total running time  = " << avg_time << " us" << std::endl;
}

int main(int argc, char *argv[])
{
    unsigned loop_times = 10;
    unsigned num_vertices = 0;
    if(argc >= 4)
    {
        const char* test_mode = argv[2];

        ColdEdge input = ReadInputFromTextFile(argv[1], num_vertices);

        if(std::strcmp(test_mode, "total") == 0)
        {
            auto const start_time = std::chrono::steady_clock::now();

            for (int i = 0; i < loop_times; i++)
            {
                GPU_Graph graph(num_vertices, input);
                PageRank(&graph);
                //printFinalResults(&graph);
            }  
            auto const end_time = std::chrono::steady_clock::now(); 
            PrintBenchmark(start_time, end_time, loop_times);
        }
        else if(std::strcmp(test_mode, "graph") == 0 )
        {
            auto const start_time = std::chrono::steady_clock::now();
            GPU_Graph graph(num_vertices, input);
            auto const end_time = std::chrono::steady_clock::now(); 

            PageRank(&graph);  
            PrintBenchmark(start_time, end_time, 1);          
        }
        else if(std::strcmp(test_mode, "pagerank") == 0)
        {
            GPU_Graph graph(num_vertices, input);
            auto const start_time = std::chrono::steady_clock::now();
            for (unsigned i = 0; i < loop_times; i++)
            {
                PageRank(&graph);
            }
            auto const end_time = std::chrono::steady_clock::now(); 
            PrintBenchmark(start_time, end_time, loop_times);
        }
        else
        {
            std::cout << "Invalid Input!" << std::endl;
            std::cout << "Please input the input text file name wanted in argc[1]" << std::endl;
            std::cout << "Please input the time mode(total/graph/pangerank) to be record in argc[2]" << std::endl;
            std::cout << "Please input the number of threads wanted to use in argc[3]" << std::endl;
        }
    }
    else if (argc >= 2 && argc < 4)
    {
        ColdEdge input = ReadInputFromTextFile(argv[1], num_vertices);
        GPU_Graph graph(num_vertices, input);
        auto const start_time = std::chrono::steady_clock::now();
        for (int i = 0; i < loop_times; i++)
        {
            PageRank(&graph);
            //printFinalResults(&graph);
        }  
        auto const end_time = std::chrono::steady_clock::now(); 
        PrintBenchmark(start_time, end_time, loop_times);
    }
    else
    {
        std::cout << "Invalid Input: " << std::endl;
        std::cout << "Please input the input text file name wanted in argc[1]" << std::endl;
        std::cout << "Please input the time mode(total/graph/pangerank) to be record in argc[2]" << std::endl;
        std::cout << "Please input the number of threads wanted to use in argc[3]" << std::endl;
    }
    return 0;
}
