// Copyright 2018 H-AXE
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <cctype>
#include <memory>
#include <string>
#include <utility>
#include <sstream>
#include <iomanip>

#include "glog/logging.h"

#include "base/tokenizer.h"
#include "common/engine.h"

using axe::base::Tokenizer;

using LabelType = int;
using GraphEmbeddingFormat = std::string;

template <typename T>
std::string label_to_hex_str(T i)
{
    std::stringstream stream;
    stream << std::setfill('0') << std::setw(sizeof(T) * 2)
           << std::hex << i;
    return stream.str();
}

class Edge
{
public:
    Edge() {}
    Edge(int id, int from_id, LabelType from_label, int to_id,
         LabelType to_label, LabelType label) : id_(id), from_id_(from_id), from_label_(from_label),
                                                to_id_(to_id), to_label_(to_label), label_(label) {}

    int GetId() const { return id_; }
    LabelType GetLabel() const { return label_ }
    int GetFromId() const { return from_id_; }
    int GetToId() const { return to_id_; }
    LabelType GetFromLabel() const { return from_label_; }
    LabelType GetToLabel() const { return to_label_; }

    void SetId(int id) { id_ = id; }
    void SetLabel(LabelType label) { label_ = label; }
    void SetFromId(int from_id) { from_ = from_id; }
    void SetToId(int to_id) { to_id_ = to_id; }
    void SetFromLabel(LabelType from_label) { from_label_ = from_label; }
    void SetToLabel(LabelType to_label) { to_label_ = to_label; }

    bool operator<(const Edge &other) const
    {
        if (label_ < other.label_)
            return true;
        else if (label_ == other.label_)
        {
            if (from_label_ < other.from_label_)
                return true;
            else if (from_label_ == other.from_label_)
                return to_label_ < other.to_label_;
            else
                return false;
        }
        else
            return false;
    }
    bool operator==(const Edge &other) const
    {
        return label_ == other.label_ && from_label_ == other.from_label_ && to_label_ == other.to_label_;
    }

    friend void operator<<(axe::base::BinStream &bin_stream, const Edge &e) { bin_stream << e.id_ << e.from_id_ << e.to_id_ << e.label_; }
    friend void operator>>(axe::base::BinStream &bin_stream, Edge &e) { bin_stream >> e.id_ >> e.from_id_ >> e.to_id_ >> e.label_; }

private:
    int id_;               // identification, almost useless.
    int from_id_;          // start vertex id in directed graph.
    int to_id_;            // end vertex id in directed graph.
    LabelType from_label_; // start vertex label in directed graph.
    LabelType to_label_;   // end vertex label in directed graph.
    LabelType label_;      // edge label, usually a integer value.
}

class Neighbor
{
public:
    Neighbor() {}
    Neighbor(int e_id, int v_id, LabelType v_label, LabelType e_label) : e_id_(e_id), v_id_(v_id), v_label_(v_label), e_label_(e_label) {}

    // int GetEId() const { return e_id_; }
    int GetVId() const { return v_id_; }
    LabelType GetVLabel() const { return v_label_; }
    LabelType GetELabel() const { return e_label_ }

    void SetVId(int v_id) { v_id_ = v_id; }
    void SetVLabel(LabelType v_label) { v_label_ = v_label; }
    void SetELabel(LabelType e_label) { e_label_ = e_label; }

    bool operator<(const Neighbor &other) const { return v_label_ < other.v_label_ || (v_label_ == other.v_label_ && e_label_ < other.e_label_); }
    bool operator==(const Neighbor &other) const { return v_label_ == other.v_label_ && e_label_ == other.e_label_; }

    friend void operator<<(axe::base::BinStream &bin_stream, const Neighbor &n) { bin_stream << n.e_id_ << n.v_id_ << n.v_label_ << n.e_label_; }
    friend void operator>>(axe::base::BinStream &bin_stream, Neighbor &n) { bin_stream >> n.e_id_ >> n.v_id_ >> n.v_label_ >> n.e_label_; }

private:
    int e_id_;          // identification, almost useless.
    int v_id_;          // vertex id of the neighbor.
    LabelType v_label_; // vertex label of the neighbor.
    LabelType e_label_; // label of the edge to the neighbor.
};

class Vertex
{
public:
    Vertex() : neighbors_(std::make_shared<std::vector<Neighbor>>()) {}
    Vertex(int id, LabelType label,
           const std::shared_ptr<std::vector<Neighbor>> &neighbors) : id_(id), label_(label), neighbors_(neighbors) {}

    int GetId() const { return id_; }
    LabelType GetLabel() const { return label_; }
    const std::shared_ptr<std::vector<Neighbor>> &GetNeighbors() const { return neighbors_; }
    int GetNeighborCount() const { return neighbors_->size(); }
    bool HasNeighbor(const Neighbor &neighbor) const { return std::find(neighbors_->begin(), neighbors_->end(), neighbor) != neighbors_->end(); }

    void SetId(int id) { id_ = id; }
    void SetLabel(const LabelType &label) { label_ = label; }
    void AddNeighbor(const Neighbor &neighbor) { neighbors_->push_back(neighbor); }
    void PopNeighbor() { neighbors_->pop_back(); }

    double GetMemory() const
    {
        double ret = sizeof(int) + neighbors_->size() * sizeof(int);
        return ret;
    }
    bool operator<(const Vertex &other) const { return id_ < other.id_; }
    bool operator==(const Vertex &other) const { return id_ == other.id_; }

    friend void operator<<(axe::base::BinStream &bin_stream, const Vertex &v) { bin_stream << v.id_ << v.label_ << *(v.neighbors_); }
    friend void operator>>(axe::base::BinStream &bin_stream, Vertex &v) { bin_stream >> v.id_ >> v.label_ >> *(v.neighbors_); }

private:
    int id_;                                           // identification of the vertex, must be unique in source graph.
    LabelType label_;                                  // label of the vertex.
    std::shared_ptr<std::vector<Neighbor>> neighbors_; // the adjacent list of neighbors, recording the label of neighbors and edges.
};

class Graph
{
public:
    Graph() : vertices_(std::make_shared<std::vector<Vertex>>()) {}
    Graph(int id, int support,
          const std::shared_ptr<std::vector<Vertex>> &vertices) : id_(id), support_(support), vertices_(vertices) {}

    int GetId() const { return id_; }
    int GetSupport() const { return support_; }
    bool HasVertex(const Vertex &vertex) const {return std::find(vertices_->begin(), vertices_->end(), vertex) != vertices_->end();}

    void SetId(int id) { id_ = id; }
    void SetSupport(int support) { support_ = support; }
    const std::shared_ptr<std::vector<Vertex>> &GetVertices() const { return vertices_; }
    void AddVertex(const Vertex &vertex) { vertices_->push_back(vertex); }
    
    GraphEmbeddingFormat GetGraphEmbedding() { return embedding_; }
    void EmbedGraph(int embedding_method = 0) {} // do nothing now.

    bool operator<(const Graph &other) const { return GraphEmbedding() < other.GraphEmbedding(); }
    bool operator==(const Graph &other) const { return GraphEmbedding() == other.GraphEmbedding(); }

private:
    int id_;                                        // identification, almost useless.
    int support_;                                   // the number of occurances of this graph, determine wheter subgraph is frequent or not.
    std::shared_ptr<std::vector<Vertex>> vertices_; // the vertices list of this graph.
    GraphEmbeddingFormat embedding_;                // the embedding of the graph, could be string, or someting else.
    // std::shared_ptr<std::vector<Edge>> edges_;      // the edges list of this graph, maybe useless.
};

int ReadInt(const std::string &line, size_t &ptr)
{
    int ret = 0;
    while (ptr < line.size() && !isdigit(line.at(ptr)))
        ++ptr;
    CHECK(ptr < line.size()) << "Invalid Input";
    while (ptr < line.size() && isdigit(line.at(ptr)))
    {
        ret = ret * 10 + line.at(ptr) - '0';
        ++ptr;
    }
    return ret;
};

LabelType ReadLabel(const std::string &line, size_t &ptr)
{
    return LabelType(ReadInt(line, ptr));
}

DatasetPartition<Graph> ParseLineGraph(const std::string &line)
{
    Graph g;
    size_t ptr = 0;
    Vertex vertex;
    int v_id = ReadInt(line, ptr);
    LabelType v_label = ReadLabel(line, ptr);
    vertex.SetId(v_id);
    vertex.SetLabel(v_label);
    int neighbor_count = ReadInt(line, ptr);
    for (int i = 0; i < neighbor_count; i++)
    {
        int e_id = 0;
        int e_to_id = ReadInt(line, ptr);
        LabelType e_label = ReadLabel(line, ptr);
        Neighbor neighbor(e_id, e_to_id, -1, e_label); // We don't know the label of vertex, so just give -1.
        vertex.AddNeighbor(neighbor);
    }
    g.AddVertex(vertex);
    DatasetPartition<Graph> ret;
    ret.push_back(g);
    return ret;
};

DatasetPartition<Vertex> ParseLineVertex(const std::string &line)
{
    size_t ptr = 0;
    Vertex vertex;
    int v_id = ReadInt(line, ptr);
    LabelType v_label = ReadLabel(line, ptr);
    vertex.SetId(v_id);
    vertex.SetLabel(v_label);
    int neighbor_count = ReadInt(line, ptr);
    for (int i = 0; i < neighbor_count; i++)
    {
        int e_id = 0;
        int e_to_id = ReadInt(line, ptr);
        LabelType e_label = ReadLabel(line, ptr);
        Neighbor neighbor(e_id, e_to_id, -1, e_label); // We don't know the label of vertex, so just give -1.
        vertex.AddNeighbor(neighbor);
    }
    DatasetPartition<Vertex> ret;
    ret.push_back(vertex);
    return ret;
};

using EdgeSupportPair = std::pair<Edge, int>;
using VertexSupportPair = std::pair<Vertex, int>;
using GraphSupportPair = std::pair<Graph, int>;
using GraphDataset = axe::common::Dataset<Graph>;
using GraphPartition = axe::common::DatasetPartition<Graph>;

class FSMGramiNaive : public Job
{
public:
    void Run(TaskGraph *tg, const std::shared_ptr<Properties> &config) const override
    {
        auto input = config->GetOrSet("graph", "dataset/citeseer.lg");
        int n_partitions = std::stoi(config->GetOrSet("parallelism", "1"));
        int minimal_support = std::stoi(config->GetOrSet("minimal_support", "5"));
        int n_iters = std::stoi(config->GetOrSet("n_iters", "5")); // the max number of edges in frequent subgraph.

        auto extend_with_freq_edge = [](const GraphPartition &subgraphs, const DatasetPartition<EdgeSupportPair> &freq_edges) {
            // This can be partitioned into multiple partitions, becaused subgraphs can be partitioned.

            GraphPartition ret;
            for (size_t i = 0; i < subgraphs.size(); ++i)
            {
                // for each subgraph, get a lot of extended graph based on it and freq_edges;
                for (auto &vertex : *(subgraphs[i].GetVertices()))
                {

                    for (const auto &edge_support : freq_edges)
                    {
                        auto edge = edge_support.first;
                        auto edge_freq = edge_support.second;
                        if (vertex.GetLabel() == edge.GetFromLabel())
                        {
                            // [TODO] Stuck here at 2020/12/20 1:16am.
                            // step1: Copy subgraphs[i]
                            // step2: Choose a neighbor for vertex,
                            //          the neigbor could be a vertex in subgraph, or be a new vertex
                            // step3: Create two Neighbor objects based on edge and neigboring vertex, 
                            // step4: Add the objects to their vertex's neighbors list.
                            continue;
                            // Neighbor new_neighbor(0, 0, edge.GetToLabel(), edge.GetLabel());
                            // if (!vertex.HasNeighbor(new_neighbor))
                            // {
                            //     // we can extend the subgraph with this edge.
                            //     vertex.AddNeighbor(new_neighbor);
                            //     auto new_vertices = std::make_shared<std::vector<Vertex>>(*(subgraphs[i].GetVertices()));
                            //     vertex.PopNeighbor();

                            //     Graph new_candidate(0, 0, new_vertices);
                            // } else{
                            //     // we can alse extend the subgraph by adding a new vertex.
                            //     vertex.AddNeighbor(new_neighbor);
                            //     auto new_vertices = std::make_shared<std::vector<Vertex>>(*(subgraphs[i].GetVertices()));
                            //     vertex.PopNeighbor();

                            //     Graph new_candidate(0, 0, new_vertices);
                            //     // add a new vertex to the
                            // }
                        }
                        else if (vertex.GetLabel() == edge.GetToLabel())
                        {
                            Neighbor new_neighbor(0, 0, edge.GetFromLabel(), edge.GetLabel());
                            if (!vertex.HasNeighbor(new_neighbor))
                            {
                                // we can extend the subgraph with this edge.
                            }
                        }
                    }
                }
            }

            auto src_vertices = *src_graph.at(0).GetVertices();
            for (const auto &subgraph : subgraphs)
            {
                // For each candidate subgraph, try to get a lot of extended graph from it.

                std::vector<int> sg_vids; // The vertex Id of all the vertices in the subgraph.
                for (const auto &sg_vertex : *subgraph.GetVertices())
                {
                    sg_vids.push_back(sg_vertex.GetId());
                }

                int sg_v_i = 0;
                for (const auto &sg_vertex : *subgraph.GetVertices())
                {
                    // We extend a subgraph based on its vertex by adding edge on vertex.
                    const auto &src_vertex = src_vertices.at(sg_vertex.GetId());
                    if (sg_vertex.GetNeighborCount() < src_vertex.GetNeighborCount())
                    {
                        // If the neighbour count of the subgraph-vertex is less than src vertex,
                        // that means we can extend a edge on subgraph-vertex, and then get a new candidate.
                        for (const auto &src_neighbor : *src_vertex.GetNeighbors())
                        {
                            bool edge_existed = false;
                            // If the edge in src_vertex is existing in subgraph, then we go on to next edge.
                            for (const auto &sg_neighbor : *sg_vertex.GetNeighbors())
                            {
                                if (sg_neighbor.GetELabel() == src_neighbor.GetELabel() && sg_neighbor.GetVId() == src_neighbor.GetVId())
                                {
                                    edge_existed = true;
                                    break;
                                }
                            }
                            if (!edge_existed)
                            {
                                // We can add this edge in src_vertex to subgraph.
                                auto new_sg_vertices = std::make_shared<std::vector<Vertex>>();
                                Graph new_graph(0, 1, new_sg_vertices);
                                (*new_graph.GetVertices()) = *subgraph.GetVertices(); // should be a deep copy.

                                (*new_graph.GetVertices()).at(sg_v_i).AddNeighbor(src_neighbor); // Add the new edge to the new_subgraph.
                                CHECK_EQ((*new_graph.GetVertices()).at(sg_v_i).GetNeighborCount() == *subgraph.GetVertices().at(sg_v_i).GetNeighborCount() + 1);

                                // If the neighbor vertex of this edge doesn't exist in the subgraph, add it.
                                bool neigher_existed = false;
                                for (size_t i = 0; i < sg_vids.size(); i++)
                                {
                                    if (sg_vids[i] == src_neighbor.GetVId())
                                    {
                                        neigher_existed = true;
                                        break;
                                    }
                                }
                                if (!neigher_existed)
                                {
                                    Vertex new_neighbor_vertex;
                                    new_neighbor_vertex.SetId(src_neighbor.GetVId());
                                    new_neighbor_vertex.SetLabel(src_vertices.at(src_neighbor.GetVId()).GetLabel());
                                    new_neighbor_vertex.AddNeighbor(Neighbor(src_neighbor.GetEId(), sg_v_i, src_neighbor.GetELabel()));

                                    new_graph.AddVertex(new_neighbor_vertex);
                                }

                                ret.push_back(new_graph)
                            }
                        }
                    }
                    sg_v_i++;
                }
            }
            return ret;
        };

        // the graph is composed by a list of vertex.
        auto graph_vertices = TextSourceDataset(input, tg, n_partitions)
                                  .FlatMap([](const std::string &line) { return ParseLineVertex(line); },
                                           [](const std::vector<double> &input) {
                                               double ret = 0;
                                               for (double x : input)
                                               {
                                                   ret += x;
                                               }
                                               return ret * 2;
                                           })
                                  .PartitionBy([](const Vertex &v) { return v.GetId(); }, 1);

        graph_vertices.UpdatePartition([](DatasetPartition<Vertex> &data) {
            for (auto &vertex : data)
            {
                for (auto &neighbor : *vertex.GetNeighbors())
                {
                    neighbor.SetVLabel(data.at(neighbor.GetVId()).GetLabel()); // Set the label of neighbor in vertex.
                }
            }
        });
        graph_vertices = graph_vertices.PartitionBy([](const Vertex &v) { return v.GetId(); }, n_partitions);
        graph_vertices.UpdatePartition([](DatasetPartition<Vertex> &data) {
            std::sort(data.begin(), data.end(), [](const Vertex &a, const Vertex &b) { return a.GetId() < b.GetId(); });
        }); // sort each partition with vertex id, which will make locating operations easier.

        auto graph_edges = graph_vertices.MapPartition(
            [](const DatasetPartition<Vertex> &data) {
                DatasetPartition<Edge> &ret;
                for (const auto &vertex : data)
                {
                    for (const auto &neighbor : *vertex.GetNeighbors())
                    {
                        ret.push_back(Edge(0, vertex.GetId(), vertex.GetLabel(), neighbor.GetVId(), neighbor.GetVLabel(), neighbor.GetELabel()));
                        // To get the neighbor's label, I modified the Neighbor structure, thus introduced some complexity when creating Neighbor.
                        // In fact we can get the label of neighbor by fetching in the data, (but makesure that the partition should be whole dataset)
                        // [TODO] If the label for neighbor is seldom used in Neighbor/Vertex, I think we can use the second solution.
                        // Or we can also change the input, which need additional pre-processing on dataset.
                    }
                }
                return ret;
            });

        auto edge_encoding_method = [](const Edge &e) {
            // std::vector<LabelType> vec_encoding{e.GetFromLabel(), e.GetToLabel(), e.GetLabel()};
            return label_to_hex_str(e.GetFromLabel()) + label_to_hex_str(e.GetToLabel()) + label_to_hex_str(e.GetLabel());
        }; // I was meant to define encoding method for ReduceBy(), but I found that we can define the operator== and operator< for Edge instead.

        // evaluate the support of each edge, and filter out the one whose support is less than the threshold.
        auto frequent_edges = graph_edges.MapPartition([](const DatasetPartition<Edge> &data) {
                                             DatasetPartition<EdgeSupportPair> &ret;
                                             ret.reserve(data.size());
                                             for (const auto &edge : data)
                                             {
                                                 Edge edge_pattern(0, 0, edge.GetFromLabel(), 0, edge.GetToLabel(), edge.GetELabel());
                                                 ret.push_back(std::make_pair(edge, 1)); // edge could be replaced with edge_pattern.
                                             }
                                             return ret;
                                         })
                                  .ReduceBy([](const EdgeSupportPair &esp) { return esp.first; }, [](EdgeSupportPair &esp1, EdgeSupportPair &esp2) { esp1.second += esp2.second; }, n_partitions)
                                  .MapPartition([minimal_support](const DatasetPartition<EdgeSupportPair> &data) {
                                      DatasetPartition<EdgeSupportPair> &ret;
                                      ret.reserve(data.size());
                                      for (const auto &esp : data)
                                      {
                                          if (esp.second >= minimal_support)
                                              ret.push_back(esp);
                                      }
                                      return ret;
                                  }); // Each partition will have a set of frequent edges.

        auto result = axe::common::Dataset<Graph>();

        // In the Grami paper, we will remove edge from SrcGraph and frequent_edges in each iteration.
        // This will incur two partition update in each outer loop, and our frequent edges and src_graph should be combined and partitioned each iteration.

        // Solution1: In the following section, we make the src_graph and frequent_edges be single partition datasets.
        // We present the src_graph as a list of vertices.
        auto src_graph = std::make_shared<axe::common::Dataset<Vertex>>(
            graph_vertices.PartitionBy([](const Vertex &v) { return 0 }, 1)
                .MapPartition([](DatasetPartition<Vertex> &data) {
                    GraphPartition ret;
                    // ret.push_back(Graph(0,1,data.data()));
                    auto vertices = std::make_shared<std::vector<Vertex>>();
                    for (const auto &v : data)
                        vertices->push_back(v);
                    ret.push_back(Graph(0, 1, vertices));
                    return ret;
                }));

        int frequent_edges_count = 100;
        frequent_edges = frequent_edges.PartitionBy([](const EdgeSupportPair &e) { return 0 }, 1);
        frequent_edges.UpdatePartition([&frequent_edges_count](DatasetPartition<EdgeSupportPair> &data) {
            frequent_edges_count = data.size();
            for (size_t i = 0; i < frequent_edges_count; ++i)
            {
                data.at(i).first.SetId(i); // Set a id to identify the frequent edges.
            }
        });
        // auto frequent_edges = std::make_shared<axe::common::Dataset<EdgeSupportPair>>(
        //     frequent_edges.PartitionBy([](const EdgeSupportPair &e) { return 0 }, 1));

        // outer loop, each iteration will be based on one frequent edge.
        for (size_t mining_step = 0; mining_step < frequent_edges_count; mining_step++)
        {
            auto frequent_edge_as_graph = frequent_edges.MapPartition([](const DatasetPartition<EdgeSupportPair> &data) {
                SupportPartition ret;
                auto &last_edge = data.back().first;
                auto &last_edge_frequency = data.back().second;
                Vertex v_from, v_to;
                v_from.SetId(0);
                v_to.SetId(1); // The id must be set to distinguish different vertices with same label.
                v_from.SetLabel(last_edge.GetFromLabel());
                v_to.SetLabel(last_edge.GetToLabel());
                v_from.AddNeighbor(Neighbor(0, 1, last_edge.GetToLabel(), last_edge.GetLabel())); // indrected graph.
                v_to.AddNeighbor(Neighbor(0, 0, last_edge.GetFromLabel(), last_edge.GetLabel())); // indrected graph.
                auto vertices = std::make_shared<std::vector<Vertex>>();
                vertices->push_back(v_from);
                vertices->push_back(v_to);
                ret.push_back(Graph(0, last_edge_frequency, vertices));
            });
            auto subgraphs = std::make_shared<GraphDataset>(frequent_edge_as_graph);

            // main loop, innner loop.
            for (int i = 0; i < n_iters; i++)
            {
                subgraphs = std::make_shared<GraphDataset>(
                    (*subgraphs).MapPartitionWith(&frequent_edges, extend_with_freq_edge));
            }

            src_graph.UpdatePartition([]); // remove edge;
            frequent_edges.UpdatePartition([](DatasetPartition<Edge> &data){
                data.pop_back(); // We pop_back to avoid updating the identify of the edge.
            }));
        }

        auto shared_graph = std::make_shared<axe::common::Dataset<Vertex>>(graph);

        auto graphs = graph.MapPartition([](const DatasetPartition<Vertex> &data) {
            DatasetPartition<Graph> ret;
            auto vertices = std::make_shared<std::vector<Vertex>>();
            for (auto &v : data)
                vertices->push_back(v);
            ret.push_back(Graph(0, 1, vertices));
            return ret;
        });
        auto shared_graphs = std::make_shared<axe::common::Dataset<Graph>>(graphs); // The big single graph G , with id=0, and support=1.

        auto start_candidates_label = graph.MapPartition([](const DatasetPartition<Vertex> &data) {
                                               DatasetPartition<std::pair<LabelType, int>> ret;
                                               for (const auto &v : data)
                                                   ret.push_back(std::make_pair(v.GetLabel(), 1));
                                               return ret;
                                           })
                                          .ReduceBy([](const std::pair<LabelType, int> v_c) { return v_c.first.GetLabel(); },
                                                    [](std::pair<LabelType, int> &agg, const std::pair<LabelType, int> &update) { agg.second += update.second; })
                                          .MapPartition([minimal_support](const DatasetPartition<std::pair<LabelType, int>> &data) {
                                              DatasetPartition<std::pair<LabelType, int>> ret;
                                              for (const auto &label_count : data)
                                              {
                                                  if (label_count.second >= minimal_support)
                                                  {
                                                      ret.push_back(label_count);
                                                  }
                                              }
                                              return ret;
                                          }); // currently each paration will have a set of frequent_vertex_labels with their support.
                                              //   .PartitionBy([](const std::pair<LabelType, int> &) { return 0; }, 1); // aggregate the partitions to on partition.

        auto get_graphs_with_vertex_label = [](const DatasetPartition<Vertex> &src_graph, const DatasetPartition<std::pair<LabelType, int>> &labels) {
            DatasetPartition<Graph> ret;
            for (const auto &v : src_graph)
            {
                for (const auto &label : labels)
                {
                    if (label.first == v.GetLabel())
                    {
                        Vertex new_vertex;
                        new_vertex.SetId(v.GetId());
                        new_vertex.SetLabel(v.GetLabel());

                        auto vertices = std::make_shared<std::vector<Vertex>>();
                        vertices->push_back(new_vertex);
                        ret.push_back(Graph(0, label.second, vertices));
                        break;
                    }
                }
            }
            return ret;
        };

        auto all_candidates_labels = start_candidates_label.Broadcast([](const std::pair<LabelType, int> &c) { return c.first; }, n_partitions);
        auto start_candidates = graph.SharedDataMapPartitionWith(&all_candidates_labels, get_graphs_with_vertex_label);
        auto results = start_candidates; // should be a deep copy

        auto extend_subgraph = [](const DatasetPartition<Graph> &subgraphs, const DatasetPartition<Graph> &src_graph) {
            DatasetPartition<Graph> extended_subgraphs;
            auto src_vertices = *src_graph.at(0).GetVertices();
            for (const auto &subgraph : subgraphs)
            {
                // For each candidate subgraph, try to get a lot of extended graph from it.

                std::vector<int> sg_vids; // The vertex Id of all the vertices in the subgraph.
                for (const auto &sg_vertex : *subgraph.GetVertices())
                {
                    sg_vids.push_back(sg_vertex.GetId());
                }

                int sg_v_i = 0;
                for (const auto &sg_vertex : *subgraph.GetVertices())
                {
                    // We extend a subgraph based on its vertex by adding edge on vertex.
                    // If the vertex in subgraph doesn't exist in src_graph, then we do not extend.
                    int src_graph_ptr = 0;
                    for (size_t src_graph_ptr = 0; src_graph_ptr < src_vertices.size(); ++src_graph_ptr)
                    {
                        if (src_vertices[src_graph_ptr].GetId() == sg_vertex.GetId())
                        {
                            break;
                        }
                    }
                    if (src_graph_ptr >= src_vertices.size())
                    {
                        continue;
                    }
                    const auto &src_vertex = src_vertices.at(src_graph_ptr); // Find sg_vertex.GetId() in src_vertices;

                    if (sg_vertex.GetNeighborCount() < src_vertex.GetNeighborCount())
                    {
                        // If the neighbour count of the subgraph-vertex is less than src vertex,
                        // that means we can extend a edge on subgraph-vertex, and then get a new candidate.
                        for (const auto &src_neighbor : *src_vertex.GetNeighbors())
                        {
                            bool edge_existed = false;
                            // If the edge in src_vertex is existing in subgraph, then we go on to next edge.
                            for (const auto &sg_neighbor : *sg_vertex.GetNeighbors())
                            {
                                if (sg_neighbor.GetELabel() == src_neighbor.GetELabel() && sg_neighbor.GetVId() == src_neighbor.GetVId())
                                {
                                    edge_existed = true;
                                    break;
                                }
                            }
                            if (!edge_existed)
                            {
                                // We can add this edge in src_vertex to subgraph.
                                auto new_sg_vertices = std::make_shared<std::vector<Vertex>>();
                                Graph new_graph(0, 1, new_sg_vertices);
                                (*new_graph.GetVertices()) = *subgraph.GetVertices(); // should be a deep copy.

                                (*new_graph.GetVertices()).at(sg_v_i).AddNeighbor(src_neighbor); // Add the new edge to the new_subgraph.
                                CHECK_EQ((*new_graph.GetVertices()).at(sg_v_i).GetNeighborCount() == *subgraph.GetVertices().at(sg_v_i).GetNeighborCount() + 1);

                                // Check whether the neighbor vertex existing in src_graph.
                                int src_graph_ptr = 0;
                                for (size_t src_graph_ptr = 0; src_graph_ptr < src_vertices.size(); ++src_graph_ptr)
                                {
                                    if (src_vertices[src_graph_ptr].GetId() == src_neighbor.GetVId())
                                    {
                                        break;
                                    }
                                }
                                if (src_graph_ptr < src_vertices.size()) // the neighbor vertex existing in src_graph.
                                {
                                    // If the neighbor vertex of this edge doesn't exist in the subgraph, add it.
                                    bool neigher_existed_in_sg = false;
                                    for (size_t i = 0; i < sg_vids.size(); i++)
                                    {
                                        if (sg_vids[i] == src_neighbor.GetVId())
                                        {
                                            neigher_existed_in_sg = true;
                                            break;
                                        }
                                    }
                                    if (!neigher_existed_in_sg)
                                    {
                                        Vertex new_neighbor_vertex;
                                        new_neighbor_vertex.SetId(src_neighbor.GetVId());
                                        new_neighbor_vertex.SetLabel(src_vertices.at(src_neighbor.GetVId()).GetLabel()); // Find src_neighbor.GetVId() in src_vertices;
                                        new_neighbor_vertex.AddNeighbor(Neighbor(src_neighbor.GetEId(), sg_v_i, src_neighbor.GetELabel()));

                                        new_graph.AddVertex(new_neighbor_vertex);
                                    }
                                }
                                // Add the new extended graph into extended_subgraphs, as outputs.
                                extended_subgraphs.push_back(new_graph)
                            }
                        }
                    }
                    sg_v_i++;
                }
            }
            return extended_subgraphs;
        };
        auto frequent_subgraph_isomorphism = [](const DatasetPartition<Graph> &subgraphs, const DatasetPartition<Graph> &src_graph) {
            // subgraph: a lot of candidate subgraphs
            // src_graph: the original big graph
            // progress: iterate over each subgraph, calculate its support in src_graph
            // return: the candidate that support is over minimal_support.
            DatasetPartition<Graph> frequent_subgraphs;
            return frequent_subgraphs;
        };
        auto update_results = [](DatasetPartition<Graph> &results, const DatasetPartition<Graph> &new_candidate) {
            results.AppendPartition(new_candidate);
        };

        // main loop
        auto candidates = start_candidates;
        for (size_t i = 0; i < n_iters; i++)
        {
            candidates = candidates.MapPartitionWith(shared_graphs.get(), extend_subgraph).MapPartitionWith(shared_graphs.get(), frequent_subgraph_isomorphism);
            results.UpdatePartitionWith(&candidates, update_results);
        }

        results.ApplyRead([](const DatasetPartition<Graph> &data) {
            for (auto &graph : data)
            {
                LOG(INFO) << "frequent subgraph: (" << graph.GetVertices() << ")@support=" << graph.GetSupport();
            }
            google::FlushLogFiles(google::INFO);
        });
        axe::common::JobDriver::ReversePrintTaskGraph(*tg);
    }
};
