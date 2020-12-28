# Introduction

In this project, we implemented  an niave algorithm for frequent subgraph mining (FSM) on Ursa. After that, we tried to parallelize the naive algorithm to acceralte the mining progress.

## Naive FSM algorithm

The Niave FSM algorithm consists of two stages: Extending stage and Evaluating stage.

### Extending Stage

### Evaluating Stage

## How to Parallelize FSM algortihm

Almost all the FSM algorithm following the "extend-evaluate" two-stage pattern. The first stage will extend the current frequent subgraph to get new candidates, the second stage will calculate the support for those candidates and determin whether they are frequent or not. For naive algorithm that we implemented, both extending stage and evaluating stage are time-consuming. Therefore, we decided to parallel these two stage with ursa to acceralted mining progress.

There are two methods to parallelize the Naive algorithms: 1) Parallize the Graph G 2) Parallelize the Candidates

### Parallelize Graph G

We can divide the big Graph G into multiple partitions. Each partiton only holds a part of the Graph G. In each iteration of mining, we calculate support for each partitions and reduce the results to get the final support of the subgraphs. After that we can determine who is the frequent subgraphs.

However, such methods have two drawbacks: 1) As the Graph G is partitioned, then we may miss some subgraphs which cross two or more partitions. This will lead to false negative results. 2) As we have to calculate support on each partition, and determine frequency after reducing, we will have more calculation because in non-parallel version, if we already know the support achives minimal_support, we can return directly.

### Parallelize Candidates

In each mining iteration, we will generate a set of candidates to evaluate. we can split the candidates into multiple partitions, and evaluate the frequency of their local candidates. After that, we can also extend it on that partition, and then evaluate again. This methods requires that each partiton will have a full replica of the Graph G, which will lead to too much storage overhead. However, the results is explict and we can still return directly after we know that the frequency achives mininal_frequency. 
