package graph

import (
	"DNA-Sequence-Alignments/dna_aligner/common"
	"math"
)

// FindMaximumWeightPath finds the path with maximum total weight (sum of scores of anchors in path).
// Graph: node indices are -1 (source), 0..numAnchors-1 (anchors), numAnchors (sink).
// Returns a slice of indices of anchors (0 to numAnchors-1) that form the max weight path.
func FindMaximumWeightPath(graph map[int][]common.Edge, numAnchors int) []int {
	if numAnchors == 0 {
		// Check if graph for source->sink exists, but generally no path if no anchors.
		if _, ok := graph[-1]; !ok { // If source isn't even in graph
			return []int{}
		}
	}

	dist := make(map[int]float64)
	pred := make(map[int]int) // pred[v] = u if edge u->v is in max path

	// Initialize distances: source is -1, anchors 0..N-1, sink is N
	for i := -1; i <= numAnchors; i++ {
		dist[i] = math.Inf(-1) // Negative infinity for longest path
	}
	dist[-1] = 0 // Distance from source to source is 0

	// Topological order: source, then anchors (0 to numAnchors-1), then sink.
	// Assumes anchors are processed in an order compatible with graph edges (e.g., by QueryStart).
	// The graph construction ensures this.
	topoOrder := make([]int, 0, numAnchors+2)
	topoOrder = append(topoOrder, -1) // Source
	for i := 0; i < numAnchors; i++ { // Anchor nodes
		topoOrder = append(topoOrder, i)
	}
	topoOrder = append(topoOrder, numAnchors) // Sink

	// Compute longest path using dynamic programming
	for _, u := range topoOrder {
		if uEdges, ok := graph[u]; ok { // Check if node u has outgoing edges
			for _, edge := range uEdges {
				v := edge.To
				weight := edge.Weight // Weight is the score of node v itself (as per BuildSegmentGraph)

				if dist[u] != math.Inf(-1) && dist[u]+weight > dist[v] {
					dist[v] = dist[u] + weight
					pred[v] = u
				}
			}
		}
	}

	// Reconstruct path by backtracking from sink
	var pathIndices []int
	curr := numAnchors // Start from sink

	for {
		predecessor, exists := pred[curr]
		if !exists || curr == -1 { // Reached source or no predecessor for curr
			break
		}
		// Only add actual anchor indices (0 to numAnchors-1) to the path.
		// `predecessor` is the node *before* `curr`.
		// If `curr` is an anchor node, then `predecessor` could be source or another anchor.
		// Python code logic: path.append(pred[curr]) then reverse. pred[curr] is 'u'.
		// So if edge is u->v, and curr=v, pred[v]=u. Add 'u' to path.
		// We want the nodes *in* the path, not their predecessors.
		// The path items are the anchor indices themselves.
		// If path is S -> A -> B -> T, pred[T]=B, pred[B]=A, pred[A]=S.
		// When curr=T, pred[T]=B. We add B. curr=B.
		// When curr=B, pred[B]=A. We add A. curr=A.
		// When curr=A, pred[A]=S. S is -1. We don't add S.
		// So, path is [B, A]. Reverse gives [A, B]. This is correct.
		if predecessor != -1 { // Don't add source node index to the list of anchor indices
			pathIndices = append(pathIndices, predecessor)
		}
		curr = predecessor // Move to the predecessor
	}

	// Reverse pathIndices to get correct order (start to end)
	for i, j := 0, len(pathIndices)-1; i < j; i, j = i+1, j-1 {
		pathIndices[i], pathIndices[j] = pathIndices[j], pathIndices[i]
	}

	// Filter out any potential non-anchor indices (e.g. if sink was somehow added)
	// This step should not be strictly necessary if logic is correct, but good for safety.
	finalPath := make([]int, 0, len(pathIndices))
	for _, idx := range pathIndices {
		if idx >= 0 && idx < numAnchors { // Ensure it's a valid anchor index
			finalPath = append(finalPath, idx)
		}
	}

	return finalPath
}
