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

using namespace optPagerank::baseline_graph;

extern const double damping_factor = 0.85;
extern const unsigned max_iterations = 100;
extern const double tolerance = 1e-10;

// Read Input (pairs of source and destination links) from file with format:
// src_index dest_index
// ...
// src_index dest_index
std::vector<Edge> ReadInputFromTextFile(const char *input_file, unsigned &num_vertices)
{
    std::ifstream myfile(input_file);
    std::vector<Edge> edges;
    if (myfile.is_open())
    {
        Edge e;
        while (myfile >> e.src >> e.dest)
        {
            unsigned larger = (e.src > e.dest) ? e.src : e.dest;
            num_vertices = (num_vertices > larger) ? num_vertices : larger;
            edges.push_back(e);
        }
        ++num_vertices;
        myfile.close();
    }
    return edges;
}

bool ToleranceCheck(const unsigned &num_v, Graph *graph)
{
    // Sum up the pagerank
    double pr_sum = 0.0;
    for (unsigned i = 0; i < num_v; i++)
    {
        pr_sum += graph->vertices[i]->pageRank();
    }
    // Calculate the cur_toleranceor
    pr_sum = 1.0 / pr_sum;
    double cur_tolerance = 0.0;
    for (unsigned i = 0; i < num_v; i++)
    {
        graph->vertices[i]->setPageRank(graph->vertices[i]->pageRank() * pr_sum);
        // norm 1
        cur_tolerance += std::fabs(graph->vertices[i]->pageRank() - graph->vertices[i]->prePageRank());
    }
    if (cur_tolerance < tolerance)
    {
        std::cout << "Current toleranceor: " << cur_tolerance << std::endl;
        return true;
    }
    return false;
}

void PageRank(Graph *graph)
{
    const unsigned num_v = graph->VertexesNum();
    double init_rank = double(1.0 / num_v);
    double pr_random = (1.0 - damping_factor) / num_v;

    for (unsigned i = 0; i < num_v; i++)
    {
        graph->vertices[i]->setPageRank(init_rank);
        graph->vertices[i]->setPrePageRank(0.0);
    }

    unsigned iter = 0;
    while (iter++ < max_iterations)
    {
        // Update the pagerank values in every iteration
        for (unsigned i = 0; i < num_v; i++)
        {
            graph->vertices[i]->setPrePageRank(graph->vertices[i]->pageRank());
            graph->vertices[i]->setPageRank(0.0);
        }

        // Distribute the pr_sum of all dangling nodes(no outer edges) to all nodes.
        double dangling_pr_sum = 0.0;
        for (unsigned i = 0; i < num_v; i++)
        {
            if (graph->vertices[i]->numOutwardEdges() == 0)
            {
                dangling_pr_sum += graph->vertices[i]->prePageRank();
            }
        }
        double pr_dangling = damping_factor * dangling_pr_sum / num_v;

        // Iterater all the vertexes and calculate its adjacency function l(pi,pj) of all inward links
        // Update its pagerank value by adding pr_eigenvector from its inward links separately
        for (unsigned i = 0; i < num_v; i++)
        {
            for (std::list<Vertex *>::iterator it = graph->vertices[i]->neighbors.begin();
                 it != graph->vertices[i]->neighbors.end(); ++it)
            {
                unsigned inward_edge_index = (*it)->getIndex();
                double num_outward_edges = (*it)->numOutwardEdges();
                double pr = graph->vertices[inward_edge_index]->prePageRank();

                double pr_eigenvector = damping_factor * pr / num_outward_edges;
                graph->vertices[i]->setPageRank(graph->vertices[i]->pageRank() + pr_eigenvector);
            }
            graph->vertices[i]->setPageRank(graph->vertices[i]->pageRank() + pr_random + pr_dangling);
        }

        // finish when cur_toleranceor is smaller than tolerance we set
        if (ToleranceCheck(num_v, graph))
        {
            std::cout << "Iteration time: " << iter << std::endl;
            break;
        }
    }
}

void printFinalResults(Graph *graph)
{
    std::cout << "PageRank values: \n";
    for (int i = 0; i < graph->VertexesNum(); ++i)
    {
        std::cout << "The index is: " << i << " with value " << graph->vertices[i]->pageRank() << '\n';
    }
    std::cout << '\n';
}

void PrintBenchmark(std::chrono::time_point<std::chrono::steady_clock> start_t, std::chrono::time_point<std::chrono::steady_clock> const end_t, const int loop_t)
{
    auto const avg_time = std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t).count() / double(loop_t);
    std::cout << "Average total running time  = " << avg_time << " us" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc == 3)
    {
        int loop_times = 10;
        unsigned num_vertices = 0;
        const char *test_mode = argv[2];

        std::vector<Edge> input = ReadInputFromTextFile(argv[1], num_vertices);

        if (std::strcmp(test_mode, "total") == 0)
        {
            auto const start_time = std::chrono::steady_clock::now();

            for (int i = 0; i < loop_times; i++)
            {
                Graph graph(num_vertices, input);
                PageRank(&graph);
                // printFinalResults(&graph);
            }
            auto const end_time = std::chrono::steady_clock::now();
            PrintBenchmark(start_time, end_time, loop_times);
        }
        else if (std::strcmp(test_mode, "graph") == 0)
        {
            auto const start_time = std::chrono::steady_clock::now();
            Graph graph(num_vertices, input);
            auto const end_time = std::chrono::steady_clock::now();

            PageRank(&graph);
            PrintBenchmark(start_time, end_time, 1);
        }
        else if (std::strcmp(test_mode, "pagerank") == 0)
        {
            Graph graph(num_vertices, input);
            auto const start_time = std::chrono::steady_clock::now();
            for (int i = 0; i < loop_times; i++)
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
