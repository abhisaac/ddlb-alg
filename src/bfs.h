#ifndef BFS_H
#define BFS_H

#include<iostream>
#include <list>
#include <algorithm>
#include <vector>
#include <iterator>
using namespace std;

// This class represents a directed graph using adjacency list representation
class Graph
{
    int V;    // No. of vertices
    list<int> *adj;    // Pointer to an array containing adjacency lists
public:
    Graph(int V);  // Constructor
    void addEdge(int v, int w); // function to add an edge to graph
    void BFS(int s, int Radius, vector<int> &nodes);  // prints BFS traversal from a given source s
};

#endif // BFS_H
