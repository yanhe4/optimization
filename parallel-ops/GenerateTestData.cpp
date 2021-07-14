#include <iostream>
#include <vector>   // std::vector
#include <fstream>  // std::ofstream
#include <random> 	// std::rand, std::srand, std::default_random_engine
#include <algorithm> // std::generate

// Generate a integer from min_num to max_num with uniform distribution
int numGenerator(double min_num, double max_num) 
{
   std::random_device rand_dev;
   std::mt19937 generator(rand_dev());
   std::uniform_int_distribution<int> distr(0, max_num);
   return distr(generator);
}

// Generate a integer from min_num to max_num with normal distribution
int  numNorGenerator(double min_num, double max_num) 
{
   std::random_device rand_dev;
   std::mt19937 generator(rand_dev());
   double mean = 0.5 * (max_num-min_num) + min_num;
   double stddev = 0.2 * (max_num-min_num);
   std::normal_distribution<double> distribution(mean, stddev);
   double ret;
   do
   {
      ret = distribution(generator);
   } while((ret<min_num)||(ret>=max_num));
   return int(ret);
}

int main (int argc, char *argv[])
{
	if(argc == 3)
	{
      std::ofstream inputFile;
      inputFile.open("input.txt");

      int num_edges = 0, cur_edges = 0;
      int num_vertices = atoi(argv[1]);
      num_edges = atoi(argv[2]);
      // matrix to store flags of edges(already exist/empty), so that we can avoid repeating edges
      std::vector<std::vector<char>> matrix(num_vertices, std::vector<char>(num_vertices,0));
      while(cur_edges < num_edges) 
      {
      	// x are sources which are randomly generated
      	// y are destinations, some of which are generated with normal distribution
        int x = numGenerator(0, num_vertices-1);
        int y;
        // Let 20% of links generated in a specific area
        if( cur_edges < (num_edges/5) )
        {
          y = numNorGenerator(0.0, double(num_vertices-1));
        } 
        else
        {
           y = numGenerator(0.0, double(num_vertices-1));
        }
        if((x != y) && (matrix[x][y] == 0))
        {
           matrix[x][y] = 1;
           ++cur_edges;
           inputFile << x << " " << y << '\n';
         //   std::cout<< x << " " << y <<'\n';
        }
      }
      inputFile.close();
   } 
   else 
   {
      std::cerr<< "Command line not valid, please follow: \n"
               << "programName verticesNum edgesNum\n"
               << "Example: \n"
               << "./generator 100 1000 >output.txt\n";
   }
  return 0;
}

