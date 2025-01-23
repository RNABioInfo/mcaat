## Conceptual summary of the `CycleFinder` class and its functions:

1. **CycleFinder Constructor**: Initializes the `CycleFinder` object with the given parameters and starts the cycle finding process.

2. **_IncomingNotEqualToCurrentNode**: Checks if any incoming edge of a node is not equal to the node itself.

3. **_BackgroundCheck**: Performs a background check on a neighbor node to determine if it meets certain criteria.

4. **_GetOutgoings**: Retrieves the outgoing edges of a node that pass the background check.

5. **_GetIncomings**: Retrieves the incoming edges of a node that pass the background check.

6. **RelaxLock**: Relaxes the lock values for nodes during backtracking to allow further exploration.

7. **ProcessNeighbors**: Processes the neighbors of the current node and updates the path and stack accordingly.

8. **Backtrack**: Performs backtracking during the cycle finding process to explore alternative paths.

9. **FindCycle**: Finds cycles starting from a given node using a depth-first search approach with backtracking.

10. **FindCycleUtil**: Utility function to initialize and start the cycle finding process from a given node.

11. **ChunkStartNodes**: Chunks the start nodes based on their multiplicity for parallel processing.

12. **DepthLevelSearch**: Performs a depth-limited search to determine if a path exists between two nodes within a certain depth.

13. **FindCycles**: Finds all cycles in the graph by iterating over chunked start nodes and utilizing parallel processing.