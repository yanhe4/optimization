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

using namespace CSC586C::optimize_graph;

extern const double damping_factor = 0.85;
extern const unsigned max_iterations = 100;
extern const double tolerance = 1e-10;

// Read Input (pairs of source and destination links) from file with format:
// src_index dest_index
// ... 
// src_index dest_index 
std::vector<ColdEdge> ReadInputFromTextFile(const char* input_file, unsigned& num_vertices)
{
    std::ifstream myfile (input_file);
    std::vector<ColdEdge> edges;
    if (myfile.is_open()) 
    {
      ColdEdge e;
      while(myfile >> e.src >> e.dest)
      {
         unsigned larger = (e.src > e.dest)? e.src : e.dest;
         num_vertices = (num_vertices > larger)? num_vertices : larger;
         edges.push_back(e);
      }
      ++num_vertices;
      myfile.close();
    }
    return edges;
}

bool ToleranceCheck(const unsigned& num_v, std::vector<HotData>& nodes)
{
    // Sum up the pagerank
    double pr_sum = 0.0;
    for (unsigned i = 0; i < num_v; i++) 
    {
        pr_sum += nodes[i].pagerank;
    }
    // Calculate the cur_toleranceor
    pr_sum = 1.0 / pr_sum;
    double cur_tolerance = 0.0;
    for (unsigned i = 0; i < num_v; i++)
    {
        nodes[i].pagerank *= pr_sum;
        // norm 1
        cur_tolerance += std::fabs(nodes[i].pagerank - nodes[i].pre_pagerank);
    }
    if (cur_tolerance < tolerance)
    {
        std::cout << "Current toleranceor: " << cur_tolerance << std::endl;
        return true;
    }
    return false;
}

void PageRank(Optimized_Graph *graph)
{
    const unsigned num_v = graph->VertexesNum();
    double init_rank = double(1.0 / num_v);
    double pr_random = (1.0 - damping_factor) / num_v;
    
    for (unsigned i = 0; i < num_v; i++)
    {
        graph->nodes[i].pagerank = init_rank;
        graph->nodes[i].pre_pagerank = 0.0;
    }

    unsigned iter = 0;
    while (iter++ < max_iterations)
    {
        double dangling_pr_sum = 0.0;
        // Update the pagerank values in every iteration
        for (unsigned i = 0; i < num_v; i++)
        {
            graph->nodes[i].pre_pagerank = graph->nodes[i].pagerank;
            graph->nodes[i].pagerank = 0.0;
            
            if(graph->nodes[i].outgoing_edges_num == 0) 
            {
                dangling_pr_sum += graph->nodes[i].pre_pagerank;
            }
        }

        double pr_dangling = damping_factor * dangling_pr_sum / num_v;

        // Iterater all the vertexes and calculate its adjacency function l(pi,pj) of all inward links
        // Update its pagerank value by adding pr_eigenvector from its inward links separately
        for (unsigned i = 0; i < num_v; i++)
        {
            unsigned inward_edges_num = graph->adjEdges[i].size();

            for (unsigned edge_num = 0; edge_num < inward_edges_num; edge_num++)
            {
                unsigned inward_edge_index = graph->adjEdges[i].at(edge_num);

                double pr_eigenvector = damping_factor * graph->nodes[inward_edge_index].pre_pagerank 
                                        / graph->nodes[inward_edge_index].outgoing_edges_num;
                graph->nodes[i].pagerank += pr_eigenvector;
            }
            graph->nodes[i].pagerank += (pr_random + pr_dangling);
        }

        // finish when cur_toleranceor is smaller than tolerance we set
        if(ToleranceCheck(num_v, graph->nodes)) 
        {
            std::cout << "Iteration time: " << iter << std::endl;
            break;
        }
    }
}

void printFinalResults(Optimized_Graph* graph)
{
    std::cout << "PageRank values: \n";
    for(int i = 0; i < graph->VertexesNum(); ++i)
    {
        std::cout << "The index is: " << i << " with value " << graph->nodes[i].pagerank << '\n';  
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
    if(argc == 3)
    {
        unsigned loop_times = 10;
        unsigned num_vertices = 0;
        const char* test_mode = argv[2];

        std::vector<ColdEdge> input = ReadInputFromTextFile(argv[1], num_vertices);

        if(std::strcmp(test_mode, "total") == 0)
        {
            auto const start_time = std::chrono::steady_clock::now();

            for (int i = 0; i < loop_times; i++)
            {
                Optimized_Graph graph(num_vertices, input);
                PageRank(&graph);
                // printFinalResults(&graph);
            } 

            auto const end_time = std::chrono::steady_clock::now(); 
            PrintBenchmark(start_time, end_time, loop_times);
        }
        else if(std::strcmp(test_mode, "graph") == 0 )
        {
            auto const start_time = std::chrono::steady_clock::now();
            Optimized_Graph graph(num_vertices, input);
            auto const end_time = std::chrono::steady_clock::now(); 

            PageRank(&graph);  
            PrintBenchmark(start_time, end_time, 1);          
        }
        else if(std::strcmp(test_mode, "pagerank") == 0)
        {
            Optimized_Graph graph(num_vertices, input);
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
        }
    }
    else
    {
        std::cout << "Invalid Input: " << std::endl;
        std::cout << "Please input the input text file name wanted in argc[1]" << std::endl;
        std::cout << "Please input the time mode(total/graph/pangerank) to be record in argc[2]" << std::endl;
    }

    return 0;
}
