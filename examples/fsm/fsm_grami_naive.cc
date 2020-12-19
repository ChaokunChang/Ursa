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

class Vertex
{
public:
    Vertex() : neighbors_(std::make_shared<std::vector<int>>()) {}
    Vertex(int id, std::string label,
           const std::shared_ptr<std::vector<int>> &neighbors) : id_(id), label_(label), neighbors_(neighbors) {}

    int GetId() const { return id_; }
    std::string GetLabel() const { return label_; }
    const std::shared_ptr<std::vector<int>> &GetNeighbors() const { return neighbors_; }
    int GetNeighborCount() const { return adj_.size(); }

    void SetId(int id) { id_ = id; }
    void SetLabel(const std::string &label) { label_ = label; }
    void AddNeighbor(int neighbor_id) { neighbors_.push_back(neighbor_id); }

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
    std::string label_;
    std::shared_ptr<std::vector<int>> neighbors_; // store the vetex.id of neighbors.
};

class Edge
{
public:
    Edge() {}
    Edge(int id, std
         : string label, int from, int to) : id_(id), label_(label), from_(from), to_(to) {}

    int GetId() const { return id_; }
    int GetLabel() const { return label_ }
    int GetFrom() const { return from_; }
    int GetTo() const { return to_; }

    bool operator<(const Edge &other) const { return id_ < other.id_; }
    bool operator==(const Edge &other) const { return id_ == other.id_; }

    friend void operator<<(axe::base::BinStream &bin_stream, const Edge &e) { bin_stream << e.id_ << e.label_ << e.from_ << e.to_; }
    friend void operator>>(axe::base::BinStream &bin_stream, Edge &e) { bin_stream >> e.id_ >> e.label_ >> e.from_ >> e.to_; }

private:
    int id_;
    int from_;
    int to_;
    std::string label_;
}

class Graph
{
public:
    Vertex() : vertices_(std::make_shared<std::vector<Vertex>>()), edges_(std::make_shared<std::vector<Edge>>()) {}
    Vertex(int id, int support,
           const std::shared_ptr<std::vector<Vertex>> &vertices,
           const std::shared_ptr<std::vector<Edge>> &edges) : id_(id), support_(support), vertices_(vertices), edges_(edges) {}

    int GetId() const { return id_; }
    int GetSupport() const { return support_; }
    const std::shared_ptr<std::vector<Vertex>> &GetVertices() const { return vertices_; }
    const std::shared_ptr<std::vector<Edge>> &GetEdges() const { return edges_; }
    void AddVertex(const Vertex &v) { vertices_.push_back(v); }
    void AddEdge(const Edge &e) { edges_.push_back(e) }
    bool Empty() const { return vertices_.empty() && edges_.empty(); }

private:
    int id_;
    int support_;
    std::shared_ptr<std::vector<Vertex>> vertices_;
    std::shared_ptr<std::vector<Edge>> edges_;
}

int
ReadInt(const std::string &line, size_t &ptr)
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
}

void ParseLine(const std::string &line, Graph &g)
{
    size_t ptr = 0;
    if (line[0] == '#')
    {
        // New graph. Create a new Graph object.
        // g = new Graph(); do nothing.
    }
    else if (line[0] == 'v')
    {
        // New vertex.
        ptr++;
        int v_id = ReadInt(line, ptr);
        int v_label = ReadInt(line, ptr);
        Vertex new_vertex;
        new_vertex.SetId(v_id);
        new_vertex.SetLabel(v_label);
        g.AddVertex(new_vertex);
    }
    else if (line[0] == 'e')
    {
        // New edge.
        ptr++;
        int e_id = ReadInt(line, ptr);
        int e_from = ReadInt(line, ptr);
        int e_to = ReadInt(line, ptr);
        int e_label = ReadInt(line, ptr);
        Edge new_edge(e_id, e_label, e_from, e_to);
        g.AddEdge(new_edge);
        auto vertices = g.GetVertices();
        CHECK(vertices.size() > e_from && vertices.size() > e_to) << "invalid edge";
        vertices->at(e_from).AddNeighbor(e_to);
        vertices->at(e_to).AddNeighbor(e_from);
    }
    else
    {
        CHECK(false) << "Invalid input."
    }
}


