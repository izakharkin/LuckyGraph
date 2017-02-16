//
// Created by ilya on 16.02.17.
// Copyright (c) 2017, ilya. All rights reserved.
//

#ifndef LUCKYGRAPH_VERTEX_H
#define LUCKYGRAPH_VERTEX_H


class Vertex {
 public:
  int index_;
  Vertex()
      : index_(-1) {}
  explicit Vertex(int index)
      : index_(index) {}
};

class WeightedVertex : public Vertex {
 public:
  double weight_;

  WeightedVertex()
      : Vertex(), weight_(-1) {}
  WeightedVertex(int index, double weight)
      : Vertex(index), weight_(weight) {}
};


#endif //LUCKYGRAPH_VERTEX_H
