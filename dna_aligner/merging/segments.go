package merging

import (
	"DNA-Sequence-Alignments/dna_aligner/common"
	"DNA-Sequence-Alignments/dna_aligner/config"
	"math"
)

// MergeAdjacentSegments merges adjacent or nearly adjacent segments.
// Input segments MUST be sorted by QueryStart.
// Segments are (q_start, q_end, r_start, r_end) inclusive.
func MergeAdjacentSegments(segments []common.Segment, maxGapUser int) []common.Segment {
	if len(segments) <= 1 { // No merging needed for 0 or 1 segment
		// Return a copy to avoid modifying input if it's a slice header from elsewhere
		result := make([]common.Segment, len(segments))
		copy(result, segments)
		return result
	}

	// Ensure segments are sorted (caller's responsibility, but good to be aware)
	// For safety, if this function can be called with unsorted segments:
	// sortedSegments := make([]common.Segment, len(segments))
	// copy(sortedSegments, segments)
	// sort.Slice(sortedSegments, func(i, j int) bool { ... })
	// For now, assume `segments` is already sorted.

	merged := []common.Segment{segments[0]} // Start with the first segment

	for i := 1; i < len(segments); i++ {
		currentMerged := &merged[len(merged)-1] // Pointer to the last segment in `merged`
		nextSegToConsider := segments[i]

		// Calculate gaps (qGap can be negative if there's overlap)
		qGap := nextSegToConsider.QueryStart - currentMerged.QueryEnd - 1
		rGap := nextSegToConsider.RefStart - currentMerged.RefEnd - 1

		// Python's merge condition:
		// (q_gap <= max_gap and r_gap <= max_gap and
		//  abs(q_gap - r_gap) <= max(5, min(q_gap, r_gap) * MAX_GAP_RATIO_DIFFERENCE))

		canMerge := false
		if qGap <= maxGapUser && rGap <= maxGapUser {
			minActualGapVal := math.Min(float64(qGap), float64(rGap))
			// maxDiffAllowed: max of 5 or a ratio of the smaller gap.
			// If minActualGapVal is negative (overlap), its product with ratio is also negative.
			// max(5, negative_value) = 5. So for overlaps, diff must be <=5.
			maxDiffAllowed := math.Max(5.0, minActualGapVal*config.MaxGapRatioDifference)
			if math.Abs(float64(qGap-rGap)) <= maxDiffAllowed {
				canMerge = true
			}
		}

		if canMerge {
			// Merge: extend currentMerged segment's end to nextSegToConsider's end
			currentMerged.QueryEnd = nextSegToConsider.QueryEnd
			currentMerged.RefEnd = nextSegToConsider.RefEnd
		} else {
			// No merge, add nextSegToConsider as a new segment to the merged list
			merged = append(merged, nextSegToConsider)
		}
	}
	return merged
}
