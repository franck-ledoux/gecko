/*----------------------------------------------------------------------------*/
#ifndef GMDS_MCTS_GRAPH_H
#define GMDS_MCTS_GRAPH_H
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/utils/Exception.h>
/*----------------------------------------------------------------------------*/
#include <limits.h>
#include <set>
#include <map>
#include <memory>
/*----------------------------------------------------------------------------*/
namespace gecko {
/*----------------------------------------------------------------------------*/
namespace cblock {
/*----------------------------------------------------------------------------*/
// A structure to represent a
// node in adjacency list
struct  AdjListNode
{
  gmds::TCellID dest;
  double weight;
  std::shared_ptr<AdjListNode> next;
  AdjListNode(gmds::TCellID ANodeID=0, double AWeight=0, std::shared_ptr<AdjListNode> AN=nullptr);
};
/*----------------------------------------------------------------------------*/
// A structure to represent
// an adjacency list
struct  AdjList
{
  // Pointer to head node of list
  std::shared_ptr<AdjListNode> head;
  AdjList( std::shared_ptr<AdjListNode> h=nullptr);
};
/*----------------------------------------------------------------------------*/
// A structure to represent a graph.
// A graph is an array of adjacency lists.
// Size of array will be V (number of
// vertices in graph)
class  Graph
{
public:
  Graph(const std::set<gmds::TCellID>& ANodeIDs);
  //Be careful, input ids are in the global input space numbering, not the local one
  void addEdge(gmds::TCellID ASrcNode, gmds::TCellID ADestNode, double AWeight);

  // The main function that calculates
  // distances of shortest paths from src to all
  // vertices. It is a O(ELogV) function
  void computeDijkstra(gmds::TCellID ASrcNode);
  //Be careful, input ids are in the global input space numbering, not the local one
 	void setWeight(const gmds::TCellID AN1, const gmds::TCellID AN2, const double AW);
	/**
	 *
	 * @return for each node of the graph D, you get the path from the source node S to D.
	 * The first item of the vector is D, and the last one S.
	 */
  std::map<gmds::TCellID , std::vector<gmds::TCellID> > getShortestPath();
  std::map<gmds::TCellID , double > getShortestPathWeights();
private:
  std::vector<std::pair<gmds::TCellID , double> > getAdjNodes(const gmds::TCellID  ANodeId);
  void buildShortestPaths(const gmds::TCellID ASrc);
  // A utility function to create
  // a new adjacency list node
  void  newAdjListNode(gmds::TCellID ADest, double AWeight, std::shared_ptr<AdjList> ANode);
  void  setWeightOneWay(const gmds::TCellID ASrc, const gmds::TCellID ADest, const double AWeight);

private:
  int m_nb_vertices;
  std::set<gmds::TCellID> m_input_vertices;
  std::vector<gmds::TCellID> m_graph_to_input_vertices;
  std::map<gmds::TCellID , gmds::TCellID> m_input_to_graph_vertices;
  std::vector<std::shared_ptr<AdjList>> m_array;

  std::vector<double> m_dist;

  std::map<gmds::TCellID , std::vector<gmds::TCellID> > m_shortest_path_nodes;
  std::map<gmds::TCellID , double > m_shortest_path_weight;
};
/*----------------------------------------------------------------------------*/
}     // namespace mctsblock
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_MCTS_GRAPH_H


