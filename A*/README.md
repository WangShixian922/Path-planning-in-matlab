The A* algorithm is easy to understand provided that breadth-first search (BFS) is studied. The difference between them is the rule of popping the node from the 
queue (A* calls it an open set). To be more specific, A* pops the node with the minimum f(n) while BFS follows the FIFO.

The main process of A* is maintaining the open and closed set. For each search loop, all the f values of the current node’s successors are updated, which will be put in
to the open set. After searching, the current node is marked as closed. Once the target point is included in successors, the searching process is done. And the optimal
path can be obtained by a backing search of the predecessor, starting from the target point.
