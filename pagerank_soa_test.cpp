#include <iostream>
#include <fstream> // std::ifstream
#include <sstream> // std::stringstream
#include <string>  // std::string, std::stoi
#include <cstring> // std::strcmp
#include <cmath>
#include <vector>
#include <chrono>
#include <ctime>
#include "Graph.hpp"

using namespace optPagerank::soa_graph;

extern const double damping_factor = 0.85;
extern const unsigned max_iterations = 100;
extern const double tolerance = 1e-10;

// Read Input (pairs of source and destination links) from file with format:
// src_index dest_index
// ...
// src_index dest_index
ColdEdge ReadInputFromTextFile(const char *input_file, unsigned &num_vertices)
{
    std::ifstream myfile(input_file);
    ColdEdge edges;
    unsigned source, destination;
    if (myfile.is_open())
    {
        while (myfile >> source >> destination)
        {
            unsigned larger = (source > destination) ? source : destination;
            num_vertices = (num_vertices > larger) ? num_vertices : larger;
            edges.src.push_back(source);
            edges.dest.push_back(destination);
        }
        ++num_vertices;
        myfile.close();
    }
    return edges;
}

bool ToleranceCheck(const unsigned &num_v, HotData &hotData)
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
        // std::cout << "Current toleranceor: " << cur_tolerance << std::endl;
        return true;
    }
    return false;
}

void PageRank(SoA_Graph *graph)
{
    const unsigned num_v = graph->VertexesNum();
    double init_rank = double(1.0 / num_v);
    double pr_random = (1.0 - damping_factor) / num_v;

    for (unsigned i = 0; i < num_v; i++)
    {
        graph->hotData.pagerank[i] = init_rank;
        graph->hotData.pre_pagerank[i] = 0.0;
    }

    unsigned iter = 0;
    while (iter++ < max_iterations)
    {
        double dangling_pr_sum = 0.0;
        // Update the pagerank values in every iteration
        for (unsigned i = 0; i < num_v; i++)
        {
            graph->hotData.pre_pagerank[i] = graph->hotData.pagerank[i];
            graph->hotData.pagerank[i] = 0.0;

            if (graph->hotData.outgoing_edges_num[i] == 0)
            {
                dangling_pr_sum += graph->hotData.pre_pagerank[i];
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

                double pr_eigenvector = damping_factor * graph->hotData.pre_pagerank[inward_edge_index] / graph->hotData.outgoing_edges_num[inward_edge_index];
                graph->hotData.pagerank[i] += pr_eigenvector;
            }
            graph->hotData.pagerank[i] += (pr_random + pr_dangling);
        }

        // finish when cur_toleranceor is smaller than tolerance we set
        if (ToleranceCheck(num_v, graph->hotData))
        {
            // std::cout << "Iteration time: " << iter << std::endl;
            break;
        }
    }
}

void printFinalResults(SoA_Graph *graph)
{
    std::cout << "PageRank values: \n";
    for (int i = 0; i < graph->VertexesNum(); ++i)
    {
        std::cout << "The index is: " << i << " with value " << graph->hotData.pagerank[i] << '\n';
    }
    std::cout << '\n';
}

void PrintBenchmark(std::chrono::time_point<std::chrono::steady_clock> start_t, std::chrono::time_point<std::chrono::steady_clock> const end_t, const unsigned loop_t)
{
    auto const avg_time = std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t).count() / double(loop_t);
    std::cout << "Average total running time  = " << avg_time << " us" << std::endl;
}

int main(int argc, char *argv[])
{
    int loop_times = 10;
    unsigned num_vertices = 0;
    ColdEdge input = ReadInputFromTextFile(argv[1], num_vertices);
    for (int i = 0; i < loop_times; i++)
    {
        SoA_Graph graph(num_vertices, input);
        PageRank(&graph);
        // printFinalResults(&graph);
    }
}