//
// Created by ilya on 16.02.17.
// Copyright (c) 2017, ilya. All rights reserved.
//

#include "Bfs.h"

void BFS::Run(int start_vertex) {
  pstat_->Init();
  colors_.resize(pgraph_->GetNodeCount(), ColorType::kWhite);
  bfs(start_vertex);
}

void BFS::bfs(int start_vertex) {
  queue<int> state_queue;
  state_queue.push(start_vertex);
  while (!state_queue.empty()) {
    int cur_vertex = state_queue.front();
    state_queue.pop();
    for (auto neighbour : pgraph_->GetNeighbours(cur_vertex)) {
      if (colors_[neighbour.index_] == ColorType::kWhite) {
        state_queue.push(neighbour.index_);
        colors_[neighbour.index_] = ColorType::kGrey;
        parent_[neighbour.index_] = cur_vertex;
      }
    }
    colors_[cur_vertex] = ColorType::kBlack;
  }
}