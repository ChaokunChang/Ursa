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

#include "glog/logging.h"

#include "base/tokenizer.h"
#include "common/engine.h"

using axe::base::Tokenizer;

using LabelType = int;

class Neighbor
{
public:
    Neighbor() {}
    Neighbor(int e_id, int v_id, LabelType e_label) : v_id_(v_id), e_id_(e_id), e_label_(e_label) {}

    int GetVId() const { return v_id_; }
    int GetEId() const { return e_id_; }
    LabelType GetELabel() const { return e_label_ }

    bool operator<(const Neighbor &other) const { return e_id_ < other.e_id_; }
    bool operator==(const Neighbor &other) const { return e_id_ == other.e_id_; }

    friend void operator<<(axe::base::BinStream &bin_stream, const Neighbor &n) { bin_stream << n.e_id_ << n.v_id_ << n.e_label_; }
    friend void operator>>(axe::base::BinStream &bin_stream, Neighbor &n) { bin_stream >> n.e_id_ >> n.v_id_ >> n.e_label_; }

private:
    int e_id_;
    int v_id_;
    LabelType e_label_;
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

    void SetId(int id) { id_ = id; }
    void SetLabel(const LabelType &label) { label_ = label; }
    void AddNeighbor(const Neighbor &neighbor) { neighbors_->push_back(neighbor); }

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
    int id_;
    LabelType label_;
    std::shared_ptr<std::vector<Neighbor>> neighbors_; // store the vetex.id of neighbors.
};

class Graph
{
public:
    Graph() : vertices_(std::make_shared<std::vector<Vertex>>()) {}
    Graph(int id, int support,
          const std::shared_ptr<std::vector<Vertex>> &vertices) : id_(id), support_(support), vertices_(vertices) {}

    int GetId() const { return id_; }
    int GetSupport() const { return support_; }
    void SetId(int id) { id_ = id; }
    void SetSupport(int support) { support_ = support; }
    const std::shared_ptr<std::vector<Vertex>> &GetVertices() const { return vertices_; }
    void AddVertex(const Vertex &vertex) { vertices_->push_back(vertex); }

private:
    int id_;
    int support_;
    std::shared_ptr<std::vector<Vertex>> vertices_;
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
        int e_to = ReadInt(line, ptr);
        LabelType e_label = ReadLabel(line, ptr);
        Neighbor neighbor(e_id, e_to, e_label);
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
        int e_to = ReadInt(line, ptr);
        LabelType e_label = ReadLabel(line, ptr);
        Neighbor neighbor(e_id, e_to, e_label);
        vertex.AddNeighbor(neighbor);
    }
    DatasetPartition<Vertex> ret;
    ret.push_back(vertex);
    return ret;
};

class FSMNaiveParallel : public Job
{
public:
    void Run(TaskGraph *tg, const std::shared_ptr<Properties> &config) const override
    {
        auto input = config->GetOrSet("graph", "dataset/citeseer.lg");
        int n_partitions = std::stoi(config->GetOrSet("parallelism", "5"));
        int minimal_support = std::stoi(config->GetOrSet("minimal_support", "5"));
        int n_iters = std::stoi(config->GetOrSet("n_iters", "5")); // the max number of edges in frequent subgraph.
        // means the graph is composed by a list of vertex.
        auto graph = TextSourceDataset(input, tg, n_partitions)
                         .FlatMap([](const std::string &line) { return ParseLineVertex(line); },
                                  [](const std::vector<double> &input) {
                                      double ret = 0;
                                      for (double x : input)
                                      {
                                          ret += x;
                                      }
                                      return ret * 2;
                                  })
                         .PartitionBy([](const Vertex &v) { return v.GetId(); }, n_partitions);
        graph.UpdatePartition([](DatasetPartition<Vertex> &data) {
            std::sort(data.begin(), data.end(), [](const Vertex &a, const Vertex &b) { return a.GetId() < b.GetId(); });
        }); // sort each partition with vertex id, which will make locating operations easier.
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