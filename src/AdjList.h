//
// Created by ilya on 16.02.17.
// Copyright (c) 2017, ilya. All rights reserved.
//

#ifndef LUCKYGRAPH_ADJLIST_H
#define LUCKYGRAPH_ADJLIST_H

#include "Vertex.h"

class AdjList {
 public:
  AdjList();
  AdjList(const vector<vector<int>>& adj_matrix);
  size_t Size() const;
  bool HasNegativeEdges() const;
  vector<WeightedVertex> GetNeighbours(size_t vertex_num) const;
 private:
  vector<vector<WeightedVertex>> adj_list_;
  bool has_negative_edges_;
};

#endif //LUCKYGRAPH_ADJLIST_H
