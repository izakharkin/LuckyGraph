//
// Created by ilya on 16.02.17.
// Copyright (c) 2017, ilya. All rights reserved.
//

#include "Statistics.h"

EdgesStat::EdgesStat()
    : has_negative_edges_(false) {}

void EdgesStat::Init() {}

void EdgesStat::NegativeEdgeDetected() {
  has_negative_edges_ = true;
}

bool EdgesStat::HasNegativeEdges() const {
  return has_negative_edges_;
}


ComponentStat::ComponentStat()
    : node_stack_(), stack_size_(0), components_(0) {}

void ComponentStat::Init() {}

void ComponentStat::EnterNode() {
  if (node_stack_.empty()) {
    components_ += 1;
    stack_size_ += 1;
  }
}

int ComponentStat::GetComponentsCount() const {
  return components_;
}


CycleStat::CycleStat()
    : cycle_type_(CycleType::kAcyclic), cycle_sign_(CycleSign::kNone) {}

void CycleStat::Init() {}

bool CycleStat::IsCyclic() const {
  return cycle_type_ == CycleType::kCyclic;
}

bool CycleStat::CycleFound(const CycleSign& cycle_sign) {
  cycle_sign_ = cycle_sign;
}