//
// Created by ilya on 16.02.17.
// Copyright (c) 2017, ilya. All rights reserved.
//

#include "EdgeList.h"

EdgeList::EdgeList()
    : edge_list_(0) {}

EdgeList::EdgeList(const vector<vector<int>>& adj_matrix) {
  size_t num_of_vertices = adj_matrix.size();
  edge_list_.resize(num_of_vertices);
  for (int i = 0; i < num_of_vertices; ++i) {
    for (int j = 0; j < num_of_vertices; ++j) {
      if (i != j) {
        double weight = adj_matrix[i][j];
        edge_list_.push_back(Edge(Vertex(i), Vertex(j), weight));
      }
    }
  }
}

EdgeList::EdgeList(const AdjList& adj_list) {
  for (int i = 0; i < adj_list.Size(); ++i) {
    for (auto item : adj_list.GetNeighbours(i)) {
      edge_list_.push_back(Edge(Vertex(i), Vertex(item.index_), item.weight_));
    }
  }
}

size_t EdgeList::Size() const {
  return edge_list_.size();
}

vector<Edge> EdgeList::GetEdges() const {
  return edge_list_;
}
