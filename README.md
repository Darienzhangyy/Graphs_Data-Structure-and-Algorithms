[![wercker status](https://app.wercker.com/status/3c63a22240cddc9211c6e4a8d771f374/m "wercker status")](https://app.wercker.com/project/bykey/3c63a22240cddc9211c6e4a8d771f374)

http://www2.stat.duke.edu/~cr173/Sta523_Fa15/hw/hw1.html

Function - is_valid

Description - Validate the graph object to ensure that it meets all requirements - Check that object is a list of lists. Check if there are names for the primary list that they are all unique. Check that each secondary list contains only edges and weights vectors that are of the appropriate type. Check that there are not any edges to non-existent vertices. Check that all weights are not less than or equal to 0. Check that every edge has a weight.

Function - is_undirected

Description - Check if the graph object is undirected, this is true if all directed edges have a complementary directed edge with the same weight in the opposite direction.

Function - is_isomorphic

Description - Check if the graph objects are isomorphic, meaning all vertices, edges, and weights are identical. Comparison of vertices should be based on names not indexes.

Function - is_connected

Description - Determine if there is any path between vertex v1 and vertex v2 in graph g. If v1 or v2 are not in g then throw an error.

Function - shortest_path

Description - Find the shortest path from vertex v1 to vertex v2 using the edges of graph g. Note that there may not be a unique solution for any given graph, you are only required to return one path.


Function - min_span_tree

Description - A tree is an undirected graph in which any two vertices are connected by exactly one path (no simple cycles). Therefore, a minimum spanning tree is the tree that connects all vertices in a graph with the shortest possible total of edges, using the existing edges. If given a directed graph return an error. Note that there may not be a unique solution for any given graph, you are only required to return one tree.

Function - plot_graph

Description - This function should be able to take any graph object and produce a reasonably attractive visual representation of that graph. Your algorithm should make use edge weights to layout the distance between vertices.
