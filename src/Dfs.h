//
// Created by ilya on 16.02.17.
// Copyright (c) 2017, ilya. All rights reserved.
//

#ifndef LUCKYGRAPH_DFS_H
#define LUCKYGRAPH_DFS_H

#include "Graph.h"

enum class ColorType {
  kWhite = 0,
  kGrey = 1,
  kBlack = 2,
};

class TimeItem {
 public:
  int time_in;
  int time_out;
};

class DFS: public ITraverse {
 public:
  DFS(unique_ptr<Graph> pgraph, unique_ptr<IStat> pstat);
  virtual void Run(int start_vertex) override;
 private:
  unique_ptr<Graph> pgraph_;
  unique_ptr<IStat> pstat_;

  vector<TimeItem> times_;
  int timer;
  vector<ColorType> colors_;
  void dfs(int cur_vertex);
};

#endif //LUCKYGRAPH_DFS_H
