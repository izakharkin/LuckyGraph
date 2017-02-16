//
// Created by ilya on 16.02.17.
// Copyright (c) 2017, ilya. All rights reserved.
//

#ifndef LUCKYGRAPH_ITRAVERSE_H
#define LUCKYGRAPH_ITRAVERSE_H

class ITraverse {
 public:
  virtual void Run(int start_vertex) = 0;
};

#endif //LUCKYGRAPH_ITRAVERSE_H
