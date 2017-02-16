//
// Created by ilya on 16.02.17.
// Copyright (c) 2017, ilya. All rights reserved.
//

#ifndef LUCKYGRAPH_STATISTICS_H
#define LUCKYGRAPH_STATISTICS_H

class StatNode {
 public:
  int data_;
};

class IStat {
 public:
  virtual void Init() = 0;
  virtual ~IStat() {}
};

class EdgesStat : public IStat {
 public:
  EdgesStat();
  void Init() override;
  void NegativeEdgeDetected();
  bool HasNegativeEdges() const;
 private:
  bool has_negative_edges_;
};

class ComponentStat : public IStat {
 public:
  ComponentStat();
  void Init() override;
  void EnterNode();
  int GetComponentsCount() const;
 private:
  stack<StatNode> node_stack_;
  int stack_size_;
  int components_;
};


enum class CycleType {
  kCyclic = 0,
  kAcyclic = 1,
};

enum class CycleSign {
  kNone = -1,
  kNegativeCycle = 0,
  kPositiveCycle = 1,
};

class CycleStat : IStat {
 public:
  CycleStat();
  void Init() override;
  bool CycleFound(const CycleSign&);
  bool IsCyclic() const;
  CycleSign GetCycleSign();
 private:
  CycleType cycle_type_;
  CycleSign cycle_sign_;
};


#endif //LUCKYGRAPH_STATISTICS_H
