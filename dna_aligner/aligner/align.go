package aligner

import (
	"DNA-Sequence-Alignments/dna_aligner/common"
	"DNA-Sequence-Alignments/dna_aligner/config"
	"DNA-Sequence-Alignments/dna_aligner/graph"
	"DNA-Sequence-Alignments/dna_aligner/matching"
	"DNA-Sequence-Alignments/dna_aligner/merging"
	"DNA-Sequence-Alignments/dna_aligner/regions"
	"DNA-Sequence-Alignments/dna_aligner/sequence"
	"fmt"
	"math"
	"sort"
)

// FindAlignment is the main alignment function.
// Returns a slice of Segments (q_start, q_end, r_start, r_end) inclusive and sorted.
func FindAlignment(query, ref string, minMatchLenUser int) []common.Segment {
	queryLen := len(query)
	refLen := len(ref)
	if queryLen == 0 || refLen == 0 {
		return []common.Segment{}
	}

	// --- Adaptive parameter selection (from Python logic) ---
	var kValuesToTry []int
	currentMinMatchLength := minMatchLenUser
	if currentMinMatchLength == 0 {
		currentMinMatchLength = config.MinMatchLength
	}
	currentStride := config.DefaultStride       // Base default, may be overridden
	currentMaxErrors := config.DefaultMaxErrors // Base default

	gcContent := sequence.CalculateGCContent(query)
	fmt.Printf("GC content: %.4f\n", gcContent)
	seqLengthConsidered := int(math.Min(float64(queryLen), float64(refLen)))

	if seqLengthConsidered < config.VeryShortSeqThreshold {
		kValuesToTry = config.VeryShortSeqKValues
		currentMinMatchLength = config.VeryShortSeqMinMatchLength
		currentStride = config.VeryShortSeqStride
		currentMaxErrors = config.HighGCMaxErrors
	} else { // Not "very short"
		if seqLengthConsidered < config.ShortSeqThreshold { // "short"
			if gcContent < config.LowGCThreshold {
				kValuesToTry = config.LowGCKValues
				currentMaxErrors = config.LowGCMaxErrors
			} else if gcContent < config.HighGCThreshold {
				kValuesToTry = config.MedGCKValues
				currentMaxErrors = config.LowGCMaxErrors + 1
			} else {
				kValuesToTry = config.HighGCKValues
				currentMaxErrors = config.HighGCMaxErrors
			}
		} else { // "long"
			if gcContent < config.LowGCThreshold {
				kValuesToTry = config.LowGCKValues
				currentMaxErrors = config.LowGCMaxErrors
			} else if gcContent < config.HighGCThreshold {
				kValuesToTry = config.MedGCKValues
				currentMaxErrors = config.LowGCMaxErrors + 1
			} else {
				kValuesToTry = config.HighGCKValues
				currentMaxErrors = config.HighGCMaxErrors
			}
		}
		// Stride for non-"very short" is k-dependent, handled in loop below.
	}
	// --- End adaptive parameters ---

	var forwardAnchors, reverseAnchors []common.AnchorMatch
	for _, k := range kValuesToTry {
		fmt.Printf("Finding anchors with k=%d...\n", k)
		iterStride := currentStride                              // Use stride determined by seq length class
		if seqLengthConsidered >= config.VeryShortSeqThreshold { // If not "very short", stride is k-dependent
			iterStride = int(math.Max(1, float64(k-5)))
		}

		fAnc := matching.FindAnchors(query, ref, k, currentMinMatchLength, iterStride, currentMaxErrors)
		rAnc := matching.FindReverseAnchors(query, ref, k, currentMinMatchLength, iterStride, currentMaxErrors)
		fmt.Printf("  Found %d forward anchors, %d reverse anchors with k=%d\n", len(fAnc), len(rAnc), k)
		forwardAnchors = append(forwardAnchors, fAnc...)
		reverseAnchors = append(reverseAnchors, rAnc...)
	}

	overlapThreshForFilter := config.HighQualityOverlapThreshold
	if gcContent < config.LowGCThreshold {
		overlapThreshForFilter += 0.02
	}
	forwardAnchors = matching.FilterAnchors(forwardAnchors, overlapThreshForFilter)
	reverseAnchors = matching.FilterAnchors(reverseAnchors, overlapThreshForFilter)
	fmt.Printf("After filtering: %d forward, %d reverse anchors remaining\n", len(forwardAnchors), len(reverseAnchors))

	// --- Process forward and reverse anchors using graph chaining ---
	var chainedFwdSegments, chainedRevSegments []common.Segment
	if len(forwardAnchors) > 0 {
		// FilterAnchors already sorts by QueryStart, which BuildSegmentGraph expects.
		fwdGraph := graph.BuildSegmentGraph(forwardAnchors)
		fwdPathIndices := graph.FindMaximumWeightPath(fwdGraph, len(forwardAnchors))
		for _, idx := range fwdPathIndices {
			anc := forwardAnchors[idx]
			chainedFwdSegments = append(chainedFwdSegments, common.Segment{anc.QueryStart, anc.QueryEnd, anc.RefStart, anc.RefEnd})
		}
	}
	if len(reverseAnchors) > 0 {
		revGraph := graph.BuildSegmentGraph(reverseAnchors)
		revPathIndices := graph.FindMaximumWeightPath(revGraph, len(reverseAnchors))
		for _, idx := range revPathIndices {
			anc := reverseAnchors[idx]
			chainedRevSegments = append(chainedRevSegments, common.Segment{anc.QueryStart, anc.QueryEnd, anc.RefStart, anc.RefEnd})
		}
	}

	// --- Combine, sort, and resolve initial overlaps ---
	combinedFromChaining := append(chainedFwdSegments, chainedRevSegments...)
	initialResolvedSegments := regions.ResolveOverlaps(combinedFromChaining) // Sorts and resolves by longer

	// --- Merge adjacent segments ---
	mergedAfterInitial := merging.MergeAdjacentSegments(initialResolvedSegments, config.AdjacentMergeMaxGap)

	// --- Ensure complete coverage ---
	// EnsureCompleteCoverage expects its input `initialSegments` to be somewhat processed (sorted, major overlaps resolved).
	// mergedAfterInitial should be sorted as MergeAdjacentSegments processes sorted input.
	segmentsAfterCoveragePass := regions.EnsureCompleteCoverage(query, ref, mergedAfterInitial)

	// --- Final merging and overlap resolution ---
	// Ensure sorted before final merge as EnsureCompleteCoverage might add segments unsortedly.
	sort.Slice(segmentsAfterCoveragePass, func(i, j int) bool {
		if segmentsAfterCoveragePass[i].QueryStart != segmentsAfterCoveragePass[j].QueryStart {
			return segmentsAfterCoveragePass[i].QueryStart < segmentsAfterCoveragePass[j].QueryStart
		}
		return segmentsAfterCoveragePass[i].RefStart < segmentsAfterCoveragePass[j].RefStart
	})
	finalMergedSegments := merging.MergeAdjacentSegments(segmentsAfterCoveragePass, config.FinalMergeMaxGap)
	finalOutputSegments := regions.ResolveOverlaps(finalMergedSegments) // Final cleanup of any overlaps

	// Clamp coordinates to sequence boundaries (Python's final step)
	clampedSegments := make([]common.Segment, 0, len(finalOutputSegments))
	for _, seg := range finalOutputSegments {
		qEndClamped := int(math.Min(float64(seg.QueryEnd), float64(queryLen-1)))
		rEndClamped := int(math.Min(float64(seg.RefEnd), float64(refLen-1)))
		// Ensure segment is still valid after clamping (start <= end)
		if seg.QueryStart <= qEndClamped && seg.RefStart <= rEndClamped {
			clampedSegments = append(clampedSegments, common.Segment{
				QueryStart: seg.QueryStart, QueryEnd: qEndClamped,
				RefStart: seg.RefStart, RefEnd: rEndClamped,
			})
		}
	}
	finalOutputSegments = clampedSegments

	// Sort for consistent output as per Python's implicit behavior / good practice
	sort.Slice(finalOutputSegments, func(i, j int) bool {
		if finalOutputSegments[i].QueryStart != finalOutputSegments[j].QueryStart {
			return finalOutputSegments[i].QueryStart < finalOutputSegments[j].QueryStart
		}
		// Tie-break by ref start, then query end, then ref end for full determinism
		if finalOutputSegments[i].RefStart != finalOutputSegments[j].RefStart {
			return finalOutputSegments[i].RefStart < finalOutputSegments[j].RefStart
		}
		if finalOutputSegments[i].QueryEnd != finalOutputSegments[j].QueryEnd {
			return finalOutputSegments[i].QueryEnd < finalOutputSegments[j].QueryEnd
		}
		return finalOutputSegments[i].RefEnd < finalOutputSegments[j].RefEnd
	})

	// Final coverage calculation (for console output)
	uncoveredInfo := regions.FindUncoveredRegions(queryLen, finalOutputSegments)
	totalUncoveredLen := 0
	for _, reg := range uncoveredInfo {
		totalUncoveredLen += (reg[1] - reg[0] + 1)
	}
	coveragePerc := 0.0
	if queryLen > 0 {
		coveragePerc = 100.0 * float64(queryLen-totalUncoveredLen) / float64(queryLen)
	}
	fmt.Printf("Final coverage: %.2f%% of query (%d segments)\n", coveragePerc, len(finalOutputSegments))

	return finalOutputSegments
}
