//
// Created by ilya on 16.02.17.
// Copyright (c) 2017, ilya. All rights reserved.
//

#include "Dfs.h"

DFS::DFS(unique_ptr<Graph> pgraph, unique_ptr<IStat> pstat)
    : pgraph_(move(pgraph)), pstat_(move(pstat)) {}

void DFS::Run(int start_vertex) {
  pstat_->Init();
  times_.resize(pgraph_->GetNodeCount());
  timer = 0;
  colors_.resize(pgraph_->GetNodeCount(), ColorType::kWhite);
  dfs(start_vertex);
}

void DFS::dfs(int cur_vertex) {
  colors_[cur_vertex] = ColorType::kGrey;
  times_[cur_vertex].time_in = timer++;
  for (const auto& next_vertex : pgraph_->GetNeighbours(cur_vertex)) {
    if (colors_[next_vertex.index_] == ColorType::kGrey) {
      if (!pgraph_->HasCycle()) {
        unique_ptr<CycleStat> cycle_stat;
        cycle_stat->CycleFound(CycleSign::kNone);
        pgraph_->UpdateCycleStat(move(cycle_stat));
      }
    } else if (colors_[next_vertex.index_] == ColorType::kWhite) {
      dfs(next_vertex.index_);
    }
  }
  colors_[cur_vertex] = ColorType::kBlack;
  times_[cur_vertex].time_out = timer++;
}


