//
// Created by ilya on 16.02.17.
// Copyright (c) 2017, ilya. All rights reserved.
//

#include "AdjList.h"

AdjList::AdjList()
    : adj_list_(0), has_negative_edges_(false) {}

AdjList::AdjList(const vector<vector<int>>& adj_matrix) {
  size_t num_of_vertices = adj_matrix.size();
  adj_list_.resize(num_of_vertices);
  for (int i = 0; i < num_of_vertices; ++i) {
    for (int j = 0; j < num_of_vertices; ++j) {
      if (i != j) {
        double weight = adj_matrix[i][j];
        if (weight < 0) {
          has_negative_edges_ = true;
        }
        adj_list_[i].push_back(WeightedVertex(j, weight));
      }
    }
  }
}

size_t AdjList::Size() const {
  return adj_list_.size();
}

bool AdjList::HasNegativeEdges() const {
  return has_negative_edges_;
}

vector<WeightedVertex> AdjList::GetNeighbours(size_t vertex_num) const {
  return adj_list_[vertex_num];
}
