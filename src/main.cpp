//
// Created by Ilya Zakharkin on 18.09.16.
//
#include <functional>
#include <iostream>
#include <memory>
#include <vector>
#include <stack>
#include <queue>
#include <set>

#include <cfloat>

#if __cplusplus == 201402L // C++14

using std::make_unique;

#else // C++11

template < typename T, typename... CONSTRUCTOR_ARGS >
        std::unique_ptr<T> make_unique( CONSTRUCTOR_ARGS&&... constructor_args )
        { return std::unique_ptr<T>( new T( std::forward<CONSTRUCTOR_ARGS>(constructor_args)... ) ); }

#endif // __cplusplus == 201402L

using std::priority_queue;
using std::make_shared;
using std::shared_ptr;
using std::unique_ptr;
using std::make_pair;
using std::function;
using std::istream;
using std::vector;
using std::queue;
using std::stack;
using std::move;
using std::pair;
using std::cout;
using std::cin;
using std::set;
using std::min;
using std::max;

const double kInf = DBL_MAX;

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

class AdjList {
 public:
  AdjList();
  AdjList(const vector<vector<int>>& adj_matrix);
  size_t Size() const;
  bool HasNegativeEdges() const;
  vector<WeightedVertex> GetNeighbours(size_t vertex_num) const;
 private:
  vector<vector<WeightedVertex>> adj_list_;
  bool has_negative_edges_;
};

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

class Edge {
 public:
  shared_ptr<Vertex> from_;
  shared_ptr<Vertex> to_;
  double weight_;

  Edge()
      : from_(), to_(), weight_(-kInf) {}
  Edge(const Vertex& from, const Vertex& to, double weight)
      : from_(make_shared<Vertex>(from)), to_(make_shared<Vertex>(to)), weight_(weight) {}
};

class EdgeList {
 public:
  EdgeList();
  EdgeList(const vector<vector<int>>& adj_matrix);
  EdgeList(const AdjList& adj_list);
  size_t Size() const;
  vector<Edge> GetEdges() const;
 private:
  vector<Edge> edge_list_;
};

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

EdgesStat::EdgesStat()
    : has_negative_edges_(false) {}

void EdgesStat::Init() {}

void EdgesStat::NegativeEdgeDetected() {
  has_negative_edges_ = true;
}

bool EdgesStat::HasNegativeEdges() const {
  return has_negative_edges_;
}

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

CycleStat::CycleStat()
    : cycle_type_(CycleType::kAcyclic), cycle_sign_(CycleSign::kNone) {}

void CycleStat::Init() {}

bool CycleStat::IsCyclic() const {
  return cycle_type_ == CycleType::kCyclic;
}

bool CycleStat::CycleFound(const CycleSign& cycle_sign) {
  cycle_sign_ = cycle_sign;
}

class Graph {
 public:
  Graph();
  ~Graph();
  Graph(Graph&& graph_to_move);
  Graph(const Graph& graph_to_copy) = delete;
  Graph(const AdjList& adj_list);
  Graph(const vector<vector<int>>& adj_matrix);

  size_t GetNodeCount() const;
  AdjList GetAdjList() const;
  vector<vector<double>> GetAdjMatrix() const; // TODO: #include <zmatrix.hpp>
  EdgeList GetEdgeList() const;

  vector<WeightedVertex> GetNeighbours(int vertex_index) const;

  void ApplyToAllEdges(const function<double(double)>& func);

  double GetWeight(int first_vertex, int second_vertex) const;

  bool HasNegativeEdges() const;

  bool HasCycle() const;
  bool HasNegativeCycle() const;

  bool IsConnected() const;
  bool IsRarefied() const;

  vector<int> GetOrder() const; // top_sort
  vector<shared_ptr<Vertex>> GetCycle() const;
  vector<Graph> GetComponents() const;
  int GetComponentsCount() const;

  void UpdateEdgesStat(unique_ptr<EdgesStat>);
  void UpdateCycleStat(unique_ptr<CycleStat>);
  void UpdateCompStat(unique_ptr<ComponentStat>);

  bool IsNetwork() const;
  // TODO: here must be flow algorithms
 private:
  int node_count_;
  int edge_count_;
  unique_ptr<AdjList> edges_;
  unique_ptr<EdgesStat> edges_stat_;
  unique_ptr<CycleStat> cycle_stat_;
  unique_ptr<ComponentStat> comp_stat_;
};

Graph::Graph()
    : node_count_(0),
      edge_count_(0),
      edges_(nullptr),
      cycle_stat_(nullptr),
      comp_stat_(nullptr) {}

Graph::Graph(Graph&& graph_to_move) {
  node_count_ = graph_to_move.node_count_;
  edge_count_ = graph_to_move.edge_count_;
  edges_ = move(graph_to_move.edges_);
  cycle_stat_ = move(graph_to_move.cycle_stat_);
  comp_stat_ = move(graph_to_move.comp_stat_);
}

Graph::~Graph() {}

Graph::Graph(const AdjList& adj_list)
    : node_count_(adj_list.Size()),
      edge_count_(),
      edges_(make_unique<AdjList>(adj_list)),
      cycle_stat_(nullptr),
      comp_stat_(nullptr)
{
  if (edges_->HasNegativeEdges()) {
    edges_stat_->NegativeEdgeDetected();
  }
}

Graph::Graph(const vector<vector<int>>& adj_matrix)
    : node_count_(adj_matrix.size()),
      edge_count_(),
      edges_(make_unique<AdjList>(adj_matrix)),
      cycle_stat_(nullptr),
      comp_stat_(nullptr) {}

size_t Graph::GetNodeCount() const {
  return node_count_;
}

vector<WeightedVertex> Graph::GetNeighbours(int vertex_index) const {
  vector<WeightedVertex> neighbours = edges_->GetNeighbours(vertex_index);
  return neighbours;
}

AdjList Graph::GetAdjList() const {
  return *edges_;
}

EdgeList Graph::GetEdgeList() const {
  return EdgeList();
}

vector<vector<double>> Graph::GetAdjMatrix() const {
  vector<vector<double>> adj_matrix;
  for (int i = 0; i < edges_->Size(); ++i) {
    for (const auto& item : edges_->GetNeighbours(i)) {
      adj_matrix[i][item.index_] = item.weight_;
    }
  }
  return adj_matrix;
}

void Graph::ApplyToAllEdges(const function<double(double)>& func) {
  for (int i = 0; i < edges_->Size(); ++i) {
    for (auto& item : edges_->GetNeighbours(i)) {
      item.weight_ = func(item.weight_);
    }
  }
}

bool Graph::HasNegativeEdges() const {
  return edges_->HasNegativeEdges();
}

bool Graph::HasCycle() const {
  return cycle_stat_->IsCyclic();
}

vector<double> FordBellman(const Graph&, int);

bool Graph::HasNegativeCycle() const {
  if (cycle_stat_ == nullptr) {
    FordBellman(*this, 0);
  }
  if (cycle_stat_->GetCycleSign() == CycleSign::kNone) {
    return false;
  } else if (cycle_stat_->GetCycleSign() == CycleSign::kNegativeCycle) {
    return true;
  }
}


bool Graph::IsConnected() const {
  return comp_stat_->GetComponentsCount() == 1;
}

bool Graph::IsRarefied() const {
  return node_count_ <= edge_count_ && edge_count_ <= 2 * node_count_;
}

int Graph::GetComponentsCount() const {
  return comp_stat_->GetComponentsCount();
}

void Graph::UpdateEdgesStat(unique_ptr<EdgesStat> new_stat) {
  edges_stat_ = move(new_stat);
}

void Graph::UpdateCycleStat(unique_ptr<CycleStat> new_stat) {
  cycle_stat_ = move(new_stat);
}

void Graph::UpdateCompStat(unique_ptr<ComponentStat> new_stat) {
  comp_stat_ = move(new_stat);
}

class GraphFactory {
 public:
  GraphFactory() = delete;
  ~GraphFactory() = delete;
  //  static unique_ptr<Graph> FromStream(istream& input_stream);
  static unique_ptr<Graph> FromAdjList(const AdjList& adj_list);
  static unique_ptr<Graph> FromAdjMatrix(const vector<vector<int>>& adj_matrix, bool is_weighted);
  static unique_ptr<Graph> MakeDirectedGraph();
  static unique_ptr<Graph> MakeUndirectedGraph();
  static unique_ptr<Graph> MakeK5();
  static unique_ptr<Graph> MakeK33();
};

unique_ptr<Graph> GraphFactory::FromAdjList(const AdjList& adj_list) {
  auto g = make_unique<Graph>(adj_list);
  return g;
}

unique_ptr<Graph> GraphFactory::FromAdjMatrix(const vector<vector<int>>& adj_matrix, bool is_weighted) {
  auto g = make_unique<Graph>(adj_matrix);
  return g;
}

class ITraverse {
 public:
  virtual void Run(int start_vertex) = 0;
};

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

vector<Vertex> AStar(const Graph& graph, int start_vertex, int end_vertex, double& distance) {
  vector<Vertex> path;
  return path;
}

class DijkstraItem {
 public:
  int index;
  double dist;

  DijkstraItem()
      : index(-1), dist(kInf) {}
  DijkstraItem(int _index, double _dist)
      : index(_index), dist(_dist) {}

  bool operator < (const DijkstraItem& second_item) const {
    return this->dist > second_item.dist;
  }
};

vector<double> Dijkstra(const Graph& graph, int start_vertex) {
  vector<double> distance(graph.GetNodeCount(), kInf);
  vector<bool> visited(graph.GetNodeCount(), 0);
  vector<int> parent(graph.GetNodeCount(), 0);
  priority_queue<DijkstraItem> state_queue;
  AdjList adj_list = graph.GetAdjList();
  distance[start_vertex] = 0.0;
  parent[start_vertex] = -1;
  state_queue.push(DijkstraItem(start_vertex, 0.0));
  DijkstraItem cur_vertex;
  int cnt = 0;
  while (!state_queue.empty() && cnt < graph.GetNodeCount()) {
    cur_vertex = state_queue.top();
    state_queue.pop();
    int from = cur_vertex.index;
    if (visited[from] == true) {
      continue;
    } else {
      visited[from] = true;
    }
    for (auto neighbour : adj_list.GetNeighbours(cur_vertex.index)) {
      int to = neighbour.index_;
      double weight = neighbour.weight_;
      if (distance[from] + weight < distance[to]) {
        distance[to] = distance[from] + weight;
        parent[to] = from;
        state_queue.push(DijkstraItem(to, distance[to]));
      }
    }
    cnt++;
  }
  return distance;
}

bool Relax(const Edge& edge, vector<double>& distance) {
  int from = edge.from_->index_;
  int to = edge.to_->index_;
  double weight = edge.weight_;
  bool relaxed = false;
  if (distance[to] > distance[from] + weight) {
    distance[to] = max(-kInf, distance[from] + weight);
    relaxed = true;
  }
  return relaxed;
}

vector<double> FordBellman(Graph& graph, int start_vertex) {
  vector<double> distance(graph.GetNodeCount(), kInf);
  distance[start_vertex] = 0.0;
  EdgeList edge_list = graph.GetEdgeList();
  vector<Edge> edges = edge_list.GetEdges();
  for (int i = 0; i < graph.GetNodeCount() - 1; ++i) {
    for (const auto &edge : edges) {
      Relax(edge, distance);
    }
  }
  for (const auto &edge : edges) {
    if (Relax(edge, distance) == true) {
      unique_ptr<CycleStat> cycle_stat;
      cycle_stat->CycleFound(CycleSign::kNegativeCycle);
      graph.UpdateCycleStat(move(cycle_stat));
    }
  }
  return distance;
}

vector<vector<double>> FloydWarshall(const Graph& graph) {
  vector<vector<double>> distances = graph.GetAdjMatrix();
  for (int k = 1; k < distances.size(); ++k) {
    for (int i = 1; i < distances.size(); ++i) {
      for (int j = 1; j < distances.size(); ++j) {
        distances[i][j] = min(distances[i][j], distances[i][k] + distances[k][j]);
      }
    }
  }
  return distances;
}

vector<vector<double>> Johnson(const Graph& graph) {
  vector<vector<double>> distances;
  return distances;
}

vector<Vertex> GetShortestPath(const Graph& graph, int start_vertex, int end_vertex, double& distance) {
  vector<Vertex> shortest_path = AStar(graph, start_vertex, end_vertex, distance);
  return shortest_path;
}

vector<double> GetShortestPaths(Graph& graph, int start_vertex) {
  vector<double> distances;
  if (!graph.HasNegativeEdges()) {
     distances = Dijkstra(graph, start_vertex);
  } else {
     distances = FordBellman(graph, start_vertex);
  }
  return distances;
}

vector<vector<double>> GetAllShortestPaths(const Graph& graph) {
  vector<vector<double>> distances;
  if (!graph.IsRarefied()) {
    distances = FloydWarshall(graph);
  } else {
    distances = Johnson(graph);
  }
  return distances;
}

int main() {
  vector<vector<int>> adj_matrix;
  int num_of_vertices, start_vertex, end_vertex;

  cin >> num_of_vertices >> start_vertex >> end_vertex;
  adj_matrix.resize(num_of_vertices, vector<int>(num_of_vertices, 0));
  for (int i = 0; i < num_of_vertices; ++i)
    for (int j = 0; j < num_of_vertices; ++j)
      cin >> adj_matrix[i][j];

  unique_ptr<Graph> graph = GraphFactory::FromAdjMatrix(adj_matrix, true);
  vector<double> paths = GetShortestPaths(*graph, start_vertex-1);

  cout << paths[end_vertex-1];

  return 0;
}