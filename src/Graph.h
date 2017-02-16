//
// Created by ilya on 16.02.17.
// Copyright (c) 2017, ilya. All rights reserved.
//

#ifndef LUCKYGRAPH_GRAPH_H
#define LUCKYGRAPH_GRAPH_H

#include "Vertex.h"

class Graph {
 public:
  Graph();
  ~Graph();
  Graph(Graph&& graph_to_move);
  Graph(const Graph& graph_to_copy) = delete;
  Graph(const AdjList& adj_list);
  Graph(const vector<vector<int>>& adj_matrix);

  size_t GetNodeCount() const;
  AdjList GetAdjList() const;
  vector<vector<double>> GetAdjMatrix() const; // TODO: #include <zmatrix.hpp>
  EdgeList GetEdgeList() const;

  vector<WeightedVertex> GetNeighbours(int vertex_index) const;

  void ApplyToAllEdges(const function<double(double)>& func);

  double GetWeight(int first_vertex, int second_vertex) const;

  bool HasNegativeEdges() const;

  bool HasCycle() const;
  bool HasNegativeCycle() const;

  bool IsConnected() const;
  bool IsRarefied() const;

  vector<int> GetOrder() const; // top_sort
  vector<shared_ptr<Vertex>> GetCycle() const;
  vector<Graph> GetComponents() const;
  int GetComponentsCount() const;

  void UpdateEdgesStat(unique_ptr<EdgesStat>);
  void UpdateCycleStat(unique_ptr<CycleStat>);
  void UpdateCompStat(unique_ptr<ComponentStat>);

  bool IsNetwork() const;
  // TODO: here must be flow algorithms
 private:
  int node_count_;
  int edge_count_;
  unique_ptr<AdjList> edges_;
  unique_ptr<EdgesStat> edges_stat_;
  unique_ptr<CycleStat> cycle_stat_;
  unique_ptr<ComponentStat> comp_stat_;
};

#endif //LUCKYGRAPH_GRAPH_H
