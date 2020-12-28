# Introduction

In this project, we implemented an naive algorithm for frequent subgraph mining (FSM) on Ursa. After that, we tried to parallelize the naive algorithm to accelerate the mining progress.

## Naive FSM algorithm

The Naive FSM algorithm consists of two stages: Extending stage and Evaluating stage. We will introduce the two stages in the following two sections.

### Extending Stage

**Input:**

* src_graph: The big Graph G to mine.
* subgraphs: The frequent subgraphs got from the last iteration. For the first iteration, a subgraph will be composed of single vertex whose label is frequent.

**Progress:**

For each subgraph, extend a vertex to get some bigger candidates.

**Output:**

A set of new candidates to be evaluated in next stage.

### Evaluating Stage

First, for each vertex in candidate graph, we find all of its matching vertices in graph G, i.e., both sharing the same vertex id. Then for each vertex, there are multiple matching vertices. We enumerate all possible matching graphs in graph G, where a map records the id of vertex in candidate graph to the id of matching vertex in graph G. After that, for each pair of vertices in candidate graph, if there is an edge connecting them, the two matching vertices in matching graph must have an edge connecting them and share the same edge id. If we find the number of isomorphism subgraphs in graph G of a candidate graph is more than minimal frequency, it can be regarded as frequent isomorphism subgraph. Note that we remove those duplicate isomorphism subgraphs, so that there won't have redundancy.

## How to Parallelize FSM algortihm

Almost all the FSM algorithm following the "extend-evaluate" two-stage pattern. The first stage will extend the current frequent subgraph to get new candidates, the second stage will calculate the support for those candidates and determine whether they are frequent or not. For naive algorithm that we implemented, both extending stage and evaluating stage are time-consuming. Therefore, we decided to parallelize these two stage with Ursa to accelerating mining progress.

There are two methods to parallelize the Naive algorithm: 1) Parallelize the Graph G 2) Parallelize the Candidates

### Parallelize Graph G

We can divide the big Graph G into multiple partitions. Each partiton only holds a part of the Graph G. In each iteration of mining, we calculate support for each partition and reduce the results to get the final support of the subgraphs. After that we can determine which are the frequent subgraphs.

However, such method have two drawbacks: 1) As the Graph G is partitioned, then we may miss some subgraphs which cross two or more partitions. This will lead to false negative results. 2) As we have to calculate support on each partition, and determine frequency after reducing, we will have more calculation because in non-parallel version, if we already know the support achives minimal_support, we can return directly.

### Parallelize Candidates

In each mining iteration, we will generate a set of candidates to evaluate. we can split the candidates into multiple partitions, and evaluate the frequency of their local candidates. After that, we can also extend it on that partition, and then evaluate again. This methods requires that each partiton have a full replica of the Graph G, which will lead to too much storage overhead. However, the results is explict and we can still return directly after we know that the frequency achives mininal_frequency. 

## How to run

### Configuration

We need 4 parameters in `.ini` file.

* parallelism: the parallelism number.
* n_iters: the number of iteration to run, also the biggest number of edges in the final frequent subgraphs.
* minimal_support: the threshold to determine whether a subgraph is frequent or not.
* graph: the location of the graph to mine.

For the `graph` dataset, we need it to be the following format:

``` txt
v.id, v.label, v.neighbors.size(), neighbor_0.vid, neighbor_0.elabel, ... neighbor_i.vid, neighbor_i.elabel, ...
...
```

### Results

The frequent subgraphs will be printed in log files. here is an example.

![algorithm outputs](https://github.com/ChaokunChang/Ursa/blob/master/report_fig_1.png)