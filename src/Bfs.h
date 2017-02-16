//
// Created by ilya on 16.02.17.
// Copyright (c) 2017, ilya. All rights reserved.
//

#ifndef LUCKYGRAPH_BFS_H
#define LUCKYGRAPH_BFS_H

#include "Graph.h"

class BFS: public ITraverse {
 public:
  BFS(unique_ptr<Graph> pgraph, unique_ptr<IStat> pstat)
      : pgraph_(move(pgraph)), pstat_(move(pstat)) {}
  virtual void Run(int start_vertex) override;
 private:
  unique_ptr<Graph> pgraph_;
  unique_ptr<IStat> pstat_;

  vector<ColorType> colors_;
  vector<int> parent_;
  void bfs(int start_vertex);
};

#endif //LUCKYGRAPH_BFS_H
