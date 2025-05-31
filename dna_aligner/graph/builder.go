package graph

import "DNA-Sequence-Alignments/dna_aligner/common"

// BuildSegmentGraph builds a directed acyclic graph from anchor matches.
// Input anchors MUST be sorted by QueryStart.
// Graph nodes: -1 (source), 0..len(anchors)-1 (anchor indices), len(anchors) (sink).
func BuildSegmentGraph(anchors []common.AnchorMatch) map[int][]common.Edge {
	numAnchors := len(anchors)
	graph := make(map[int][]common.Edge)

	// Source node index: -1
	// Sink node index: numAnchors

	// Edges from source to all anchor nodes
	sourceEdges := make([]common.Edge, 0, numAnchors)
	for i := 0; i < numAnchors; i++ {
		sourceEdges = append(sourceEdges, common.Edge{To: i, Weight: anchors[i].Score})
	}
	graph[-1] = sourceEdges

	// Edges between anchor nodes and from anchor nodes to sink
	for i := 0; i < numAnchors; i++ {
		anchorIqEnd := anchors[i].QueryEnd // Inclusive end from AnchorMatch
		currentAnchorEdges := make([]common.Edge, 0)

		for j := i + 1; j < numAnchors; j++ {
			anchorJqStart := anchors[j].QueryStart // Inclusive start
			// If anchorJ starts after anchorI ends (no overlap in query)
			if anchorJqStart > anchorIqEnd {
				currentAnchorEdges = append(currentAnchorEdges, common.Edge{To: j, Weight: anchors[j].Score})
			}
		}
		// Edge from anchor i to sink (weight 0, sink doesn't add to score itself)
		currentAnchorEdges = append(currentAnchorEdges, common.Edge{To: numAnchors, Weight: 0})
		graph[i] = currentAnchorEdges
	}

	// Ensure sink node exists in graph map, even if no anchors point to it (e.g., numAnchors is 0)
	if _, ok := graph[numAnchors]; !ok {
		graph[numAnchors] = []common.Edge{} // Sink has no outgoing edges
	}
	if numAnchors == 0 { // Special case: no anchors, ensure source has edge to sink
		graph[-1] = []common.Edge{{To: 0 /*sink index when numAnchors=0*/, Weight: 0}}
	}

	return graph
}
