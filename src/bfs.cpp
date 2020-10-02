// Program to print BFS traversal from a given source vertex. BFS(int s)
// traverses vertices reachable from s.
#include "bfs.h"

Graph::Graph(int V)
{
  this->V = V;
  adj = new list <int>[V];
}

void Graph::addEdge(int v, int w)
{
  adj[v].push_back(w); // Add w to v’s list.
  adj[w].push_back(v);
}

void Graph::BFS(int s, int Radius, vector<int>& nodes)
{
  // Mark all the vertices as not visited
  bool *visited = new bool[V];
  int *radius = new int[V];
  for(int i = 0; i < V; i++){
    visited[i] = false;
    radius[i] = 0;
  }

  // Create a queue for BFS
  list<int> queue;

  // Mark the current node as visited and enqueue it
  visited[s] = true;
  queue.push_back(s);

  // 'i' will be used to get all adjacent vertices of a vertex
  while(!queue.empty())
  {
    // Dequeue a vertex from queue and print it
    s = queue.front();
    queue.pop_front();

    // Get all adjacent vertices of the dequeued vertex s
    // If a adjacent has not been visited, then mark it visited
    // and enqueue it

    for(auto i : adj[s]) //i = adj[s].begin(); i != adj[s].end(); ++i)
    {

      if(!visited[i])
      {
        radius[i] = radius[s]+1;
        visited[i] = true;
        queue.push_back(i);
      }
    }
  }

  for (int i = 0; i < V; ++i) {
    if(radius[i] <= Radius)
      nodes.push_back(i);
  }
}
