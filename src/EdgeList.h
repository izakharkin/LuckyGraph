//
// Created by ilya on 16.02.17.
// Copyright (c) 2017, ilya. All rights reserved.
//

#ifndef LUCKYGRAPH_EDGELIST_H
#define LUCKYGRAPH_EDGELIST_H

#include "Vertex.h"

class Edge {
 public:
  shared_ptr<Vertex> from_;
  shared_ptr<Vertex> to_;
  double weight_;

  Edge()
      : from_(), to_(), weight_(-kInf) {}
  Edge(const Vertex& from, const Vertex& to, double weight)
      : from_(make_shared<Vertex>(from)), to_(make_shared<Vertex>(to)), weight_(weight) {}
};

class EdgeList {
 public:
  EdgeList();
  EdgeList(const vector<vector<int>>& adj_matrix);
  EdgeList(const AdjList& adj_list);
  size_t Size() const;
  vector<Edge> GetEdges() const;
 private:
  vector<Edge> edge_list_;
};

#endif //LUCKYGRAPH_EDGELIST_H
