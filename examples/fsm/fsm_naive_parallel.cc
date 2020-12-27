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
#include <stack>
#include <string>
#include <utility>
#include <vector>

#include "base/tokenizer.h"
#include "common/engine.h"
#include "glog/logging.h"

using axe::base::Tokenizer;

using LabelType = int;

class Neighbor
{
public:
    Neighbor() {}
    Neighbor(int e_id, int v_id, int e_label) : v_id_(v_id), e_id_(e_id), e_label_(e_label) {}
    Neighbor(const Neighbor &other) : v_id_(other.v_id_), e_id_(other.e_id_), e_label_(other.e_label_) {}

    int GetVId() const { return v_id_; }
    int GetEId() const { return e_id_; }
    int GetELabel() const { return e_label_; }
    std::string GetStr() const { return "(n, " + std::to_string(e_id_) + "," + std::to_string(v_id_) + "," + std::to_string(e_label_) + ")"; }

    bool operator<(const Neighbor &other) const { return e_id_ < other.e_id_; }
    bool operator==(const Neighbor &other) const { return e_id_ == other.e_id_; }

    friend void operator<<(axe::base::BinStream &bin_stream, const Neighbor &n) { bin_stream << n.e_id_ << n.v_id_ << n.e_label_; }
    friend void operator>>(axe::base::BinStream &bin_stream, Neighbor &n) { bin_stream >> n.e_id_ >> n.v_id_ >> n.e_label_; }

private:
    int e_id_;
    int v_id_;
    int e_label_;
};

class Vertex
{
public:
    Vertex() : neighbors_(std::make_shared<std::vector<Neighbor>>()) {}
    Vertex(int id, int label, const std::shared_ptr<std::vector<Neighbor>> &neighbors) : id_(id), label_(label), neighbors_(neighbors) {}
    Vertex(const Vertex &other)
    {
        id_ = other.id_;
        label_ = other.label_;
        neighbors_ = std::make_shared<std::vector<Neighbor>>(*other.neighbors_);
    }

    int GetId() const { return id_; }
    int GetLabel() const { return label_; }
    const std::shared_ptr<std::vector<Neighbor>> &GetNeighbors() const { return neighbors_; }
    int GetNeighborCount() const { return neighbors_->size(); }
    std::string GetStr(bool verbose=false) const
    {
        std::string ret = "[v, " + std::to_string(id_) + "," + std::to_string(label_);
        int neighbor_count = 0;
        for (const auto &neighbor : *neighbors_)
        {
            if (!verbose && neighbor_count >= 5){
                ret += "...... (" + neighbors_->size() + "neighbors in total)";
                break;
            }
            ret += ", " + neighbor.GetStr() + " ";
            neighbor_count++;
        }
        ret += "]";
        return ret;
    }

    void SetId(int id) { id_ = id; }
    void SetLabel(const int &label) { label_ = label; }
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
    int label_;
    std::shared_ptr<std::vector<Neighbor>> neighbors_; // store the vetex.id of neighbors.
};

class Graph
{
public:
    Graph() : vertices_(std::make_shared<std::vector<Vertex>>()) {}
    Graph(int id, int support, const std::shared_ptr<std::vector<Vertex>> &vertices) : id_(id), support_(support), vertices_(vertices) {}
    Graph(const Graph &graph)
    {
        this->id_ = graph.id_;
        this->support_ = graph.support_;
        this->vertices_ = std::make_shared<std::vector<Vertex>>(*graph.vertices_);
    }

    int GetId() const { return id_; }
    int GetSupport() const { return support_; }
    std::string GetStr(bool verbose=false) const
    {
        std::string ret = "\n {g, " + std::to_string(id_) + "," + std::to_string(support_);
        int vertex_count = 0;
        for (const auto &vertex : *vertices_)
        {
            ret += "\n ********";
            if (!verbose && vertex_count >= 10){
                ret += "...... (" + vertices_->size() + "vertices in total)";
                break;
            }
            ret += ", " + vertex.GetStr(verbose) + " ";
            vertex_count++;
        }
        ret += "}";
        return ret;
    }

    void SetId(int id) { id_ = id; }
    void SetSupport(int support) { support_ = support; }
    const std::shared_ptr<std::vector<Vertex>> &GetVertices() const { return vertices_; }
    void AddVertex(const Vertex &vertex) { vertices_->push_back(vertex); }

    bool operator<(const Graph &other) const { return id_ < other.id_; }
    bool operator==(const Graph &other) const { return id_ == other.id_; }

    friend void operator<<(axe::base::BinStream &bin_stream, const Graph &g) { bin_stream << g.id_ << g.support_ << *(g.vertices_); }
    friend void operator>>(axe::base::BinStream &bin_stream, Graph &g) { bin_stream >> g.id_ >> g.support_ >> *(g.vertices_); }

private:
    int id_;
    int support_;
    std::shared_ptr<std::vector<Vertex>> vertices_;
};

using LabelSupportPair = std::pair<LabelType, int>;
using EdgeSupportPair = std::pair<Edge, int>;
using VertexSupportPair = std::pair<Vertex, int>;
using GraphDataset = axe::common::Dataset<Graph>;
using GraphPartition = axe::common::DatasetPartition<Graph>;

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

DatasetPartition<Graph> ParseLineGraph(const std::string &line)
{
    Graph g;
    size_t ptr = 0;
    Vertex vertex;
    int v_id = ReadInt(line, ptr);
    int v_label = ReadInt(line, ptr);
    vertex.SetId(v_id);
    vertex.SetLabel(v_label);
    int neighbor_count = ReadInt(line, ptr);
    for (int i = 0; i < neighbor_count; i++)
    {
        int e_id = 0;
        int e_to = ReadInt(line, ptr);
        int e_label = ReadInt(line, ptr);
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
    int v_label = ReadInt(line, ptr);
    vertex.SetId(v_id);
    vertex.SetLabel(v_label);
    int neighbor_count = ReadInt(line, ptr);
    for (int i = 0; i < neighbor_count; i++)
    {
        int e_id = 0;
        int e_to = ReadInt(line, ptr);
        int e_label = ReadInt(line, ptr);
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
        auto input = config->GetOrSet("graph", "/home/ckchang/test1.graph");
        LOG(INFO) << "Load dataset from: " << input;
        int n_partitions = std::stoi(config->GetOrSet("parallelism", "2"));
        int minimal_support = std::stoi(config->GetOrSet("minimal_support", "2"));
        int n_iters = std::stoi(config->GetOrSet("n_iters", "1")); // the max number of edges in frequent subgraph.
        // means the graph is composed by a list of vertex.
        auto graph = TextSourceDataset(input, tg, n_partitions)
                         .FlatMap([](const std::string &line) { return ParseLineVertex(line); })
                         .PartitionBy([](const Vertex &v) { return 0; }, n_partitions);
        graph.UpdatePartition([](DatasetPartition<Vertex> &data) {
            std::sort(data.begin(), data.end(), [](const Vertex &a, const Vertex &b) { return a.GetId() < b.GetId(); });
        }); // If n_partitions=1, in fact this is not needed as the data is sorted by default.

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

        (*shared_graph).ApplyRead([](const DatasetPartition<Vertex> &src_graph) {
            for (const auto &vertex : src_graph)
            {
                LOG(INFO) << "the shared graph: " << vertex.GetStr();
            }
            google::FlushLogFiles(google::INFO);
        });

        auto frequent_vertex_labels =
            graph
                .MapPartition([](const DatasetPartition<Vertex> &data) {
                    DatasetPartition<LabelSupportPair> ret;
                    for (const auto &v : data)
                        ret.push_back(std::make_pair(v.GetLabel(), 1));
                    return ret;
                })
                .ReduceBy([](const LabelSupportPair v_c) { return v_c.first; },
                          [](LabelSupportPair &agg, const LabelSupportPair &update) { agg.second += update.second; })
                .MapPartition([minimal_support](const DatasetPartition<LabelSupportPair> &data) {
                    DatasetPartition<LabelSupportPair> ret;
                    for (const auto &lc_pair : data)
                    {
                        if (lc_pair.second >= minimal_support)
                        {
                            ret.push_back(lc_pair);
                        }
                    }
                    return ret;
                })
                .PartitionBy([](const LabelSupportPair &) { return 0; }, n_partitions);
        
        frequent_vertex_labels.ApplyRead([](const DatasetPartition<LabelSupportPair> &data) {
            LOG(INFO) << "Read Frequent V labels in this Partition: " << data.size();
            for (const auto &p : data)
            {
                LOG(INFO) << "frequent_vertex_labels: " << p.first << " , " << p.second;
                google::FlushLogFiles(google::INFO);
            }
        });

        auto get_graphs_with_vertex_label = [](const DatasetPartition<Vertex> &src_graph, const DatasetPartition<LabelSupportPair> &labels) {
            LOG(INFO) << "Frequent V labels in this Partition: " << labels.size();
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
                        LOG(INFO) << "get_graphs_with_vertex_label "
                                  << "Add new vertex " << label.first;
                        break;
                    }
                }
            }
            google::FlushLogFiles(google::INFO);
            return ret;
        };
        auto all_freq_vertex_labels = std::make_shared<axe::common::Dataset<LabelSupportPair>>(frequent_vertex_labels);
        if (n_partitions > 1)
        {
            all_freq_vertex_labels = std::make_shared<axe::common::Dataset<LabelSupportPair>>(
                frequent_vertex_labels.Broadcast([](const LabelSupportPair &p) { return 0; }, n_partitions));
        }
        auto start_candidates = (*shared_graph).SharedDataMapPartitionWith(all_freq_vertex_labels.get(), get_graphs_with_vertex_label);
        start_candidates.ApplyRead([](const DatasetPartition<Graph> &data) {
            for (const auto &g : data)
            {
                LOG(INFO) << "the start candidates: " << g.GetStr();
            }
            google::FlushLogFiles(google::INFO);
        });

        auto extend_subgraph = [](const DatasetPartition<Graph> &src_graph, const DatasetPartition<Graph> &subgraphs) {
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
                                Graph new_graph(subgraph); // Should be a deep copy.
                                new_graph.SetId(0);
                                new_graph.SetId(1);
                                (*new_graph.GetVertices()).at(sg_v_i).AddNeighbor(src_neighbor); // Add the new edge to the new_subgraph.
                                CHECK_EQ((*new_graph.GetVertices()).at(sg_v_i).GetNeighborCount(), (*subgraph.GetVertices()).at(sg_v_i).GetNeighborCount() + 1);

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

                                extended_subgraphs.push_back(new_graph);
                            }
                        }
                    }
                    sg_v_i++;
                }
            }
            return extended_subgraphs;
        };
        auto frequent_subgraph_isomorphism = [minimal_support](const DatasetPartition<Graph> &src_graph, const DatasetPartition<Graph> &subgraphs) {
            // subgraph: a lot of candidate subgraphs
            // src_graph: the original big graph
            // progress: iterate over each subgraph, calculate its support in src_graph
            // return: the candidate subgraphs that support is over minimal_support.
            DatasetPartition<Graph> frequent_subgraphs;
            for (const auto &subgraph : subgraphs)
            {
                frequent_subgraphs.push_back(subgraph);
            }
            return frequent_subgraphs;
        };

        auto update_results = [](DatasetPartition<Graph> &results, const DatasetPartition<Graph> &new_candidate) {
            // results.AppendPartition(new_candidate);
            for (const auto &candidate : new_candidate)
            {
                results.push_back(candidate);
            }
        };

        // main loop
        auto candidates = std::make_shared<axe::common::Dataset<Graph>>(start_candidates);
        auto results = (*candidates).MapPartition([](const DatasetPartition<Graph> &data) {
            DatasetPartition<Graph> ret;
            for (const auto &g : data)
            {
                ret.push_back(g);
            }
            return ret;
        });
        shared_graphs = std::make_shared<axe::common::Dataset<Graph>>(
            (shared_graphs.get())->Broadcast([](Graph &g) { return 0; }, n_partitions));
        for (size_t iter = 0; iter < n_iters; iter++)
        {
            candidates = std::make_shared<axe::common::Dataset<Graph>>((*shared_graphs).SharedDataMapPartitionWith(candidates.get(), extend_subgraph));
            candidates =
                std::make_shared<axe::common::Dataset<Graph>>((*shared_graphs).SharedDataMapPartitionWith(candidates.get(), frequent_subgraph_isomorphism));
            results.UpdatePartitionWith(candidates.get(), update_results);
            results.ApplyRead([iter](const DatasetPartition<Graph> &data) {
                for (const auto &graph : data)
                {
                    LOG(INFO) << "frequent subgraph @iter" << iter << " : " << graph.GetStr();
                }
                google::FlushLogFiles(google::INFO);
            });
        }
        results.PartitionBy([](const Graph &g) { return 0; }, n_partitions)
            .ApplyRead([iter](const DatasetPartition<Graph> &data) {
                for (const auto &graph : data)
                {
                    LOG(INFO) << "final frequent subgraphs " << " : " << graph.GetStr();
                }
                google::FlushLogFiles(google::INFO);
            });

        axe::common::JobDriver::ReversePrintTaskGraph(*tg);
    }
};

int main(int argc, char **argv)
{
    axe::common::JobDriver::Run(argc, argv, FSMNaiveParallel());
    return 0;
}
