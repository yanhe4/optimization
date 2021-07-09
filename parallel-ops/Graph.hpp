
#include <iostream>
#include <vector>
#include <list>

namespace CSC586C {

// The first edition of data structure
// Baseline_graph with AoS and use std::vector<Vertex*> to store the edges and their relationships
namespace baseline_graph{
// data structure to store graph edges
struct Edge
{
	unsigned src, dest;
};

class Vertex
{
public:
   Vertex(unsigned index_) : index(index_) 
   {
      page_rank = 0.0;
      pre_pagerank = 0.0;
      num_inward_edges = 0;
      num_outward_edges = 0;
      neighbors.resize(0);
   }
   void addNeighbor(Vertex* vertex) 
   {
      neighbors.push_back(vertex);
      ++num_inward_edges;
   }
   double pageRank() const{
      return page_rank;
   }
   double prePageRank() const{
      return pre_pagerank;
   }
   void setPageRank(double page_rank_){
      page_rank = page_rank_;
   }
   void setPrePageRank(double page_rank_){
      pre_pagerank = page_rank_;
   }
   unsigned numInwardEdges() const{
      return num_inward_edges;
   }
   unsigned numOutwardEdges() const{
      return num_outward_edges;
   }
   unsigned getIndex() const{
      return index;
   }
   void outCountIncrease() {
      ++num_outward_edges;
   }
   // For debugging only
   void printVertexInfo() const
   { 
      std::cout << "Vertex index: " << index << " outward count " << num_outward_edges
                << "page rank " << page_rank << '\n'
                << "previous page rank " << pre_pagerank <<'\n';
      std::cout << "Neighbors: ";
      for(std::list<Vertex*>::const_iterator it=neighbors.begin(); it!= neighbors.end(); ++it) 
      {
         std::cout << (*it)->index << " ";
      }
      std::cout << '\n';
   }
   std::list<Vertex*> neighbors; //list of inward neighbors (sources)
private:
   double page_rank;
   double pre_pagerank;
   unsigned num_inward_edges;
   unsigned num_outward_edges;
   unsigned index;
};

class Graph
{
public:
   Graph(unsigned num_vertices_, std::vector<Edge> edges) : num_vertices(num_vertices_){
      vertices.resize(0);
      //Initialize all the vertices on this graph
      for(unsigned i=0; i<num_vertices; ++i){
         Vertex* vertex = new Vertex(i);
         vertices.push_back(vertex);
      }
      //Add all the edges info
      for(unsigned i=0; i<edges.size(); ++i){
         vertices[edges[i].dest]->addNeighbor(vertices[edges[i].src]);
         vertices[edges[i].src]->outCountIncrease();
      }
   }
   unsigned VertexesNum() const{
      return num_vertices;
   }
   void printGraphInfo() const 
   { //For debugging only
      std::cout << "==================PRINT GRAPH=====================\n";
      std::cout << "Number of vertices: "<<num_vertices<<'\n';
      double pagerank_sum = 0.0;
      for(unsigned i=0; i<num_vertices; ++i) {
         vertices[i]->printVertexInfo();
         pagerank_sum += vertices[i]->pageRank();
      }
      std::cout << "Total pagerank sum is " << pagerank_sum << '\n';
      std::cout << "=========================END======================\n";
   }
   ~Graph() {
      for(unsigned i=0; i<num_vertices; ++i){
         Vertex* vertex  = vertices[i];
         delete vertex;
      }
      vertices.resize(0);
      num_vertices = 0;
   }
   std::vector<Vertex*> vertices; //All vertices in this graph
private:
   unsigned num_vertices;
};
}// End of namespace baseline_graph





// The second edition of data structure
// Used cold/hot data and std::vector<std::vector<int> > to store the edges and their relationships
namespace optimize_graph {
// data structure to store cold data: edges (pairs of src and dest)
struct ColdEdge 
{
	unsigned int src;
	unsigned int dest;
};
// data structure to store hot data: number of outgoing links each node and its pagerank values 
struct HotData
{
	double pagerank;
	double pre_pagerank;
	int outgoing_edges_num;
};
// class to represent a graph object
class Optimized_Graph
{
public:
	// construct a vector of vectors to represent each edge index and its inward edges (hot data)
	std::vector<std::vector<int> > adjEdges;
	// construct a vector to store hot data structure
	std::vector<HotData> nodes;

	// Graph Constructor
	Optimized_Graph(unsigned Num, std::vector<optimize_graph::ColdEdge> &input): vertex_num(Num), edges(input)
	{
		adjEdges.resize(Num);
		nodes.resize(Num);
		// add edges to the directed graph
		for (auto &edge: edges)
		{
			adjEdges[edge.dest].push_back(edge.src);
			nodes[edge.src].outgoing_edges_num++;
		}
	}
	unsigned VertexesNum() {return vertex_num;}

   //For debugging only
   void printGraph() const 
   {
      std::cout << "Number of vertices is: " << vertex_num <<'\n';
      std::cout << "Hot data: \n";
      std::cout << "Pre_pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << nodes[i].pre_pagerank << ' ';
      }
      std::cout << '\n' << "Pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << nodes[i].pagerank << ' ';
      }
      std::cout << '\n' << "OutgoingEdgeCount: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << nodes[i].outgoing_edges_num << ' ';
      }
      std::cout << '\n' << "The adj matrix:\n";
      for(unsigned i = 0; i < vertex_num; ++i){
         for(unsigned j = 0; j < adjEdges[i].size(); ++j){
            std::cout << adjEdges[i][j] << "--->" << i << " ";
         }
         std::cout<<'\n';
      }
   }
   void printHotData() const 
   {
      std::cout << "Hot data: \n";
      std::cout << "Pre_pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << nodes[i].pre_pagerank <<' ';
      }
      std::cout << '\n' << "Pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << nodes[i].pagerank << ' ';
      }
   }

private:
	unsigned int vertex_num;	
	// construct a vector to store cold data structure
	std::vector<optimize_graph::ColdEdge> const &edges;
};
} // End of namespace optimize_graph





// The third edition of data structure
// Used SoA and std::vector<std::vector<int> > to store the edges and their relationships
namespace soa_graph {

// data structure to store cold data: edges (pairs of src and dest)
struct ColdEdge 
{
	std::vector<unsigned> src;
   std::vector<unsigned> dest;
};
// data structure to store hot data: number of outgoing links each node and its pagerank values 
struct HotData
{
	std::vector<double> pagerank;
	std::vector<double> pre_pagerank;
   std::vector<int> outgoing_edges_num; 
};
// class to represent a graph object
class SoA_Graph
{
public:
	// construct a vector of vectors to represent each edge index and its inward edges (hot data)
	std::vector<std::vector<int> > adjEdges;
	// construct a data structure to store hot data arrays
	HotData hotData;

	// Graph Constructor
	SoA_Graph(unsigned Num, const ColdEdge &input): vertex_num(Num), edges(input)
	{
		adjEdges.resize(Num);
      hotData.pagerank.resize(Num);
      hotData.pre_pagerank.resize(Num);
      hotData.outgoing_edges_num.resize(Num);
		// add edges to the directed graph
      for (unsigned i = 0; i < edges.src.size(); i++)
      {
         adjEdges[edges.dest[i]].push_back(edges.src[i]);
         hotData.outgoing_edges_num[edges.src[i]]++;
      }
	}
	unsigned VertexesNum() {return vertex_num;}
   //For debugging only
   void printGraph() const 
   {
      std::cout << "Number of vertices is: " << vertex_num <<'\n';
      std::cout << "Hot data: \n";
      std::cout << "Pre_pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.pre_pagerank[i] << ' ';
      }
      std::cout << '\n' << "Pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.pagerank[i] << ' ';
      }
      std::cout << '\n' << "OutgoingEdgeCount: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.outgoing_edges_num[i] << ' ';
      }
      std::cout << '\n' << "The adj graph:\n";
      for(unsigned i = 0; i < vertex_num; ++i){
         for(unsigned j = 0; j < adjEdges[i].size(); ++j){
            std::cout << adjEdges[i][j] << "--->" << i << " ";
         }
         std::cout<<'\n';
      }
   }
   void printHotData() const 
   {
      std::cout << "Hot data: \n";
      std::cout << "Pre_pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.pre_pagerank[i] <<' ';
      }
      std::cout << '\n' << "Pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.pagerank[i] << ' ';
      }
   }

private:
	unsigned int vertex_num;	
	// construct data structure to store cold data arrays
	soa_graph::ColdEdge edges;
};
} // End of namespace soa_graph




namespace optimize_matrix {
// data structure to store cold data: edges (pairs of src and dest)
struct ColdEdge 
{
   std::vector<int> src;
   std::vector<int> dest;
};
// data structure to store hot data: number of outgoing links each node and its pagerank values 
struct HotData
{
   std::vector<double> pagerank;
   std::vector<double> pre_pagerank;
   std::vector<int> outgoing_edges_num; 
};
// class to represent a graph object
class SoA_Matrix
{
public:
   // construct a vector of vectors to represent each edge index and its inward edges (hot data)
   std::vector<std::vector<int> > adjMatrix;
   // construct a data structure to store hot data arrays
   HotData hotData;

   // Matrix Constructor
   SoA_Matrix(unsigned Num, const ColdEdge &input): vertex_num(Num) {
      adjMatrix.resize(Num, std::vector<int>(Num,0));
      hotData.pagerank.resize(Num);
      hotData.pre_pagerank.resize(Num);
      hotData.outgoing_edges_num.resize(Num);
      // add edges to the adjacency matrix
      for (unsigned i = 0; i < input.src.size(); i++)
      {
         adjMatrix[input.src[i]][input.dest[i]] = 1;
         hotData.outgoing_edges_num[input.src[i]]++;
      }
   }
   unsigned VertexesNum() {return vertex_num;}

   //For debugging only
  void printMatrix() const {
     std::cout<<"Number of vertices " << vertex_num <<'\n';
     std::cout<<"Hot data: \n";
     std::cout << "Pre_pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.pre_pagerank[i] << ' ';
      }
      std::cout << '\n' << "Pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.pagerank[i] << ' ';
      }
      std::cout << '\n' << "OutgoingEdgeCount: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.outgoing_edges_num[i] << ' ';
      }
      std::cout << '\n' << "The adj matrix:\n";
      for(unsigned i = 0; i < vertex_num; ++i){
         for(unsigned j = 0; j < vertex_num; ++j){
            std::cout << adjMatrix[i][j] << " ";
         }
         std::cout<<'\n';
      }
  }
  void printHotData() const 
   {
      std::cout << "Hot data: \n";
      std::cout << "Pre_pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.pre_pagerank[i] <<' ';
      }
      std::cout << '\n' << "Pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.pagerank[i] << ' ';
      }
   }
private:
   unsigned int vertex_num;   
};
} // End of namespace optimize_matrix



// The fourth edition of data structure for GPU computng
// Used SoA and int* adjE to store the edges and their relationships
namespace gpu_graph {
// data structure to store cold data: edges (pairs of src and dest)
struct ColdEdge 
{
   std::vector<unsigned> src;
   std::vector<unsigned> dest;
};
// data structure to store hot data: number of outgoing links each node and its pagerank values 
struct HotData
{
   std::vector<double> pagerank;
   std::vector<double> pre_pagerank;
   std::vector<int> outgoing_edges_num;
};
// class to represent a graph object
class GPU_Graph
{
public:
     // construct a vector of vectors to represent each edge index and its inward edges (hot data)
     std::vector<std::vector<int> > adjEdges;

     // Store the number of inward edges of each node.
     std::vector<int> ingoing_edges_num;
     // Store every node's begin index
     // It is likely a pointer to the index of node.
     std::vector<int> beginIndex;
     // Store the inward edges in 1-D array(1-D array version of adjEdges). 
     // The size of which is the same as the number of links.
     int* adjE; 
     // Total number of edgeds
     int num_edges;

     // construct a data structure to store hot data arrays
     HotData hotData;

     // Graph Constructor
     GPU_Graph(unsigned Num, const ColdEdge &input): vertex_num(Num), edges(input)
     {
         adjEdges.resize(Num);
         hotData.pagerank.resize(Num);
         hotData.pre_pagerank.resize(Num);
         hotData.outgoing_edges_num.resize(Num);
         ingoing_edges_num.resize(Num);
         beginIndex.resize(Num);
         // add edges to the directed graph
         for (unsigned i = 0; i < edges.src.size(); i++)
         {
            adjEdges[edges.dest[i]].push_back(edges.src[i]);
            hotData.outgoing_edges_num[edges.src[i]]++;
         }
         // Initialize adjE for GPU optimization
         num_edges = edges.src.size();
         adjE = new int[num_edges];
         int index = 0;
         for(unsigned i = 0; i < Num; ++i) {
            ingoing_edges_num[i] = adjEdges[i].size();
            beginIndex[i] = index;
            for(int j=0; j<ingoing_edges_num[i]; ++j){
               adjE[index] = adjEdges[i][j];
               ++index;
            }
         }
     }
     ~GPU_Graph() {
         delete [] adjE;
     }
     unsigned VertexesNum() {return vertex_num;}

     //Print GPU related data for debugging
     void printGpuData() const
     {
         std::cout<<"========================GPU INFO BEGIN========================\n";
         for(unsigned i = 0; i < vertex_num; ++i) {
            std::cout<<" Num of ingoing edges at " << i << " is "<< ingoing_edges_num[i]
                     <<" index begins at " << beginIndex[i] << " neighbors: \n";
            int index = beginIndex[i];
            for(int j = 0; j < ingoing_edges_num[i]; ++j){
               std::cout<<adjE[index + j]<<" ";
            }
            std::cout<<'\n';
         }
         std::cout<<"========================GPU INFO END==========================\n";
     }
   //For debugging only
   void printGraph() const 
   {
      std::cout << "Number of vertices is: " << vertex_num <<'\n';
      std::cout << "Hot data: \n";
      std::cout << "Pre_pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.pre_pagerank[i] << ' ';
      }
      std::cout << '\n' << "Pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.pagerank[i] << ' ';
      }
      std::cout << '\n' << "OutgoingEdgeCount: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.outgoing_edges_num[i] << ' ';
      }
      std::cout << '\n' << "The adj graph:\n";
      for(unsigned i = 0; i < vertex_num; ++i){
         for(unsigned j = 0; j < adjEdges[i].size(); ++j){
            std::cout << adjEdges[i][j] << "--->" << i << " ";
         }
         std::cout<<'\n';
      }
   }
   void printHotData() const 
   {
      std::cout << "Hot data: \n";
      std::cout << "Pre_pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.pre_pagerank[i] <<' ';
      }
      std::cout << '\n' << "Pagerank: ";
      for(unsigned i = 0; i < vertex_num; ++i){
         std::cout << hotData.pagerank[i] << ' ';
      }
   }

private:
   unsigned int vertex_num;   
   // construct data structure to store cold data arrays
   gpu_graph::ColdEdge edges;
};
} // End of namespace gpu_graph


} // End of namespace CSC586C
