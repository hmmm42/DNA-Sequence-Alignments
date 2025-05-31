package matching

import (
	"DNA-Sequence-Alignments/dna_aligner/common"
	"DNA-Sequence-Alignments/dna_aligner/config"
	"DNA-Sequence-Alignments/dna_aligner/sequence"
	"math"
	"sort"
)

// FindAnchors finds anchor regions between query and reference.
// K, minMatchLen, stride, maxErrors: if 0, use config defaults.
func FindAnchors(query, ref string, k, minMatchLen, stride, maxErrors int) []common.AnchorMatch {
	if k == 0 {
		k = config.DefaultK
	}
	if minMatchLen == 0 {
		minMatchLen = config.MinMatchLength
	}
	if stride == 0 {
		stride = config.DefaultStride
	}
	if maxErrors == 0 {
		maxErrors = config.DefaultMaxErrors
	}
	if k <= 0 {
		return []common.AnchorMatch{}
	}

	exactMatches := FindExactMatches(query, ref, k)
	var anchors []common.AnchorMatch
	processed := make(map[[2]int]bool) // Using [2]int as map key for (q_pos, r_pos)

	for i, em := range exactMatches {
		coordKey := [2]int{em.QueryPos, em.RefPos}
		// Python: if i % stride != 0 and (q_start, r_start) in processed: continue
		// This means: process if (it's a stride hit) OR (it's not processed yet)
		if i%stride != 0 && processed[coordKey] {
			continue
		}

		anchor := ExtendMatch(query, ref, em.QueryPos, em.RefPos, em.Length, minMatchLen, maxErrors)
		if anchor != nil {
			anchors = append(anchors, *anchor)

			// Mark key positions within the new anchor as processed (Python logic)
			// Anchor struct has inclusive QueryStart, QueryEnd.
			qS, qE := anchor.QueryStart, anchor.QueryEnd
			rS, _ := anchor.RefStart, anchor.RefEnd // rE not used for marking key in python

			matchLen := (qE - qS) + 1
			strideFactor := int(math.Max(1, float64(matchLen/10)))
			if strideFactor == 0 {
				strideFactor = 1
			}

			for j := 0; j < matchLen; j += strideFactor {
				currentQ := qS + j
				currentR := rS + j                                // Assumes diagonal progression for marking
				if currentQ < len(query) && currentR < len(ref) { // Check bounds
					processed[[2]int{currentQ, currentR}] = true
				}
			}
		}
	}
	return FilterAnchors(anchors, config.DefaultOverlapThreshold)
}

// FilterAnchors filters and prioritizes anchors based on quality and overlap.
func FilterAnchors(anchors []common.AnchorMatch, overlapThreshold float64) []common.AnchorMatch {
	if len(anchors) == 0 {
		return []common.AnchorMatch{}
	}

	sort.SliceStable(anchors, func(i, j int) bool { // Sort by score (highest first)
		return anchors[i].Score > anchors[j].Score
	})

	var filtered []common.AnchorMatch
	excludedIndices := make(map[int]bool)

	for i := 0; i < len(anchors); i++ {
		if excludedIndices[i] {
			continue
		}
		anchorI := anchors[i]
		filtered = append(filtered, anchorI)

		qStartI, qEndI := anchorI.QueryStart, anchorI.QueryEnd
		rStartI, rEndI := anchorI.RefStart, anchorI.RefEnd

		for j := i + 1; j < len(anchors); j++ {
			if excludedIndices[j] {
				continue
			}
			anchorJ := anchors[j]
			qStartJ, qEndJ := anchorJ.QueryStart, anchorJ.QueryEnd
			rStartJ, rEndJ := anchorJ.RefStart, anchorJ.RefEnd

			// Calculate query overlap (inclusive coordinates)
			qOverlapStart := int(math.Max(float64(qStartI), float64(qStartJ)))
			qOverlapEnd := int(math.Min(float64(qEndI), float64(qEndJ)))
			qOverlapLen := 0
			if qOverlapEnd >= qOverlapStart {
				qOverlapLen = qOverlapEnd - qOverlapStart + 1
			}
			qLenJ := (qEndJ - qStartJ) + 1
			qOverlapRatio := 0.0
			if qLenJ > 0 {
				qOverlapRatio = float64(qOverlapLen) / float64(qLenJ)
			}

			// Calculate reference overlap
			rOverlapStart := int(math.Max(float64(rStartI), float64(rStartJ)))
			rOverlapEnd := int(math.Min(float64(rEndI), float64(rEndJ)))
			rOverlapLen := 0
			if rOverlapEnd >= rOverlapStart {
				rOverlapLen = rOverlapEnd - rOverlapStart + 1
			}
			rLenJ := (rEndJ - rStartJ) + 1
			rOverlapRatio := 0.0
			if rLenJ > 0 {
				rOverlapRatio = float64(rOverlapLen) / float64(rLenJ)
			}

			if qOverlapRatio > overlapThreshold || rOverlapRatio > overlapThreshold {
				excludedIndices[j] = true
			}
		}
	}
	sort.SliceStable(filtered, func(i, j int) bool { // Sort by query start for next steps
		return filtered[i].QueryStart < filtered[j].QueryStart
	})
	return filtered
}

// FindReverseAnchors finds anchors between query and reverse complement of reference.
func FindReverseAnchors(query, ref string, k, minMatchLen, stride, maxErrors int) []common.AnchorMatch {
	revRef := sequence.ReverseComplement(ref)
	// FindAnchors returns anchors with coordinates relative to query and revRef.
	anchorsOnRevRef := FindAnchors(query, revRef, k, minMatchLen, stride, maxErrors)

	var reverseAnchors []common.AnchorMatch
	refOriginalLen := len(ref)
	for _, anchor := range anchorsOnRevRef {
		// anchor.RefStart, anchor.RefEnd are inclusive coordinates on revRef.
		// Convert to original ref coordinates (inclusive).
		// revRef[idx] corresponds to original_ref[refOriginalLen - 1 - idx] (complemented).
		// A segment revRef[rS .. rE] corresponds to original_ref[ (L-1-rE) .. (L-1-rS) ].
		origRStart := refOriginalLen - 1 - anchor.RefEnd
		origREnd := refOriginalLen - 1 - anchor.RefStart

		reverseAnchors = append(reverseAnchors, common.AnchorMatch{
			QueryStart: anchor.QueryStart,
			QueryEnd:   anchor.QueryEnd,
			RefStart:   origRStart,
			RefEnd:     origREnd,
			Score:      anchor.Score,
			Identity:   anchor.Identity,
			// Orientation implicitly 'r', could be set if needed by consumers
		})
	}
	// FilterAnchors (called by FindAnchors) already sorts by QueryStart.
	return reverseAnchors
}
