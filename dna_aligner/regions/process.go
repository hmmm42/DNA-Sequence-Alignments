package regions

import (
	"DNA-Sequence-Alignments/dna_aligner/common"
	"DNA-Sequence-Alignments/dna_aligner/config"
	"DNA-Sequence-Alignments/dna_aligner/matching"

	"math"
	"sort"
)

// FindUncoveredRegions finds regions in the query not covered by segments.
// Segments are (q_start, q_end, r_start, r_end) with inclusive query coords.
// Input segments should be sorted by QueryStart.
func FindUncoveredRegions(queryLen int, segments []common.Segment) [][2]int { // Returns list of [start, end] inclusive
	if len(segments) == 0 {
		if queryLen > 0 {
			return [][2]int{{0, queryLen - 1}}
		}
		return [][2]int{}
	}

	// Ensure segments are sorted by query start position (caller should ideally ensure this)
	sortedSegments := make([]common.Segment, len(segments))
	copy(sortedSegments, segments) // Avoid modifying original slice if passed by ref
	sort.Slice(sortedSegments, func(i, j int) bool {
		return sortedSegments[i].QueryStart < sortedSegments[j].QueryStart
	})

	var uncovered [][2]int
	currentPos := 0 // Marks the end of the last covered region + 1

	for _, seg := range sortedSegments {
		qStart, qEnd := seg.QueryStart, seg.QueryEnd // inclusive
		if qStart > currentPos {
			uncovered = append(uncovered, [2]int{currentPos, qStart - 1})
		}
		currentPos = int(math.Max(float64(currentPos), float64(qEnd+1)))
	}

	if currentPos < queryLen {
		uncovered = append(uncovered, [2]int{currentPos, queryLen - 1})
	}
	return uncovered
}

// findMatchesInRegionCore is a helper for finding matches in a given query region.
// Coords in returned AnchorMatch are relative to queryRegion string.
func findMatchesInRegionCore(queryRegion, ref string, minMatchL, maxErr int, kValuesToTry []int, stride int) []common.AnchorMatch {
	var regionMatches []common.AnchorMatch

	for _, k := range kValuesToTry {
		if k <= 0 || k > len(queryRegion) || k > len(ref) { // k must be valid and within bounds
			continue
		}

		// Forward anchors
		fAnchors := matching.FindAnchors(queryRegion, ref, k, minMatchL, stride, maxErr)
		for _, anc := range fAnchors {
			m := anc // Make a copy to set orientation
			m.Orientation = 'f'
			regionMatches = append(regionMatches, m)
		}

		// Reverse anchors
		rAnchors := matching.FindReverseAnchors(queryRegion, ref, k, minMatchL, stride, maxErr)
		for _, anc := range rAnchors {
			m := anc // Make a copy
			m.Orientation = 'r'
			regionMatches = append(regionMatches, m)
		}
	}
	return regionMatches
}

// FindMatchesInRegion finds matches for a smaller query region against the full reference.
// Coords in returned AnchorMatch are relative to queryRegion string.
func FindMatchesInRegion(queryRegion, ref string, minMatchLenUser, maxErrorsUser int) []common.AnchorMatch {
	segmentLen := len(queryRegion)
	if segmentLen == 0 {
		return []common.AnchorMatch{}
	}

	var kSize int
	if segmentLen < 50 {
		kSize = config.SmallKForShortSegments
	} else if segmentLen < 100 {
		kSize = config.MediumKForShortSegments
	} else {
		kSize = config.LargeKForShortSegments
	}
	kValuesToTry := []int{kSize}

	minMatchL := minMatchLenUser
	if minMatchL == 0 {
		minMatchL = config.MinMatchLength // Use global default
	}

	maxErr := maxErrorsUser
	if maxErr == 0 { // Python uses max(3, min_match_length // 10)
		maxErr = int(math.Max(3, float64(minMatchL/10)))
		if maxErr == 0 {
			maxErr = 1
		} // Ensure at least 1
	}

	return findMatchesInRegionCore(queryRegion, ref, minMatchL, maxErr, kValuesToTry, 1) // Stride 1 as in Python
}

// FindMatchesInLargeRegion finds multiple matches for a large query region using a divide-and-conquer approach.
// Coords in returned AnchorMatch are relative to queryRegion string.
func FindMatchesInLargeRegion(queryRegion, ref string, maxSegSizeUser, minMatchLenUser, maxErrorsUser int) []common.AnchorMatch {
	regionLen := len(queryRegion)
	if regionLen == 0 {
		return []common.AnchorMatch{}
	}

	maxSegSize := maxSegSizeUser
	if maxSegSize == 0 {
		maxSegSize = config.MaxSegmentSize
	}
	minMatchL := minMatchLenUser
	if minMatchL == 0 {
		minMatchL = config.MinMatchLength
	}
	maxErrDefault := maxErrorsUser
	if maxErrDefault == 0 {
		maxErrDefault = config.DefaultMaxErrors
	}

	if regionLen <= maxSegSize { // Not "large" enough, use simpler method
		return FindMatchesInRegion(queryRegion, ref, minMatchL, maxErrDefault)
	}

	var matches []common.AnchorMatch
	chunkSize, overlap := 0, 0

	if regionLen > 5000 {
		chunkSize = config.LargeRegionChunkSize
		overlap = config.LargeRegionOverlap
	} else {
		chunkSize = maxSegSize
		overlap = int(float64(chunkSize) / config.StandardChunkOverlapRatio)
	}
	if chunkSize <= 0 {
		chunkSize = config.MaxSegmentSize
	}
	if overlap < 0 {
		overlap = 0
	}
	if overlap >= chunkSize { // Ensure overlap is smaller than chunk, and positive
		overlap = chunkSize / 3
		if overlap <= 0 {
			overlap = 1
		}
	}
	if chunkSize-overlap <= 0 { // Ensure positive step
		if chunkSize > 1 {
			overlap = chunkSize - 1
		} else {
			overlap = 0
		}
	}

	for i := 0; i < regionLen; i += (chunkSize - overlap) {
		chunkStart := i
		chunkEnd := int(math.Min(float64(i+chunkSize), float64(regionLen)))

		if chunkEnd-chunkStart < minMatchL {
			continue
		}
		chunk := queryRegion[chunkStart:chunkEnd]
		chunkLen := len(chunk)

		var kValues []int
		if chunkLen < 100 {
			kValues = config.VeryShortSegmentKValues
		} else if chunkLen < 300 {
			kValues = config.ShortSegmentKValues
		} else {
			kValues = config.LongerSegmentKValues
		}

		currentMaxErrors := maxErrDefault
		if i == 0 || chunkEnd == regionLen { // Boundary chunks
			currentMaxErrors += config.BoundaryExtraErrors
		}

		chunkMatches := findMatchesInRegionCore(chunk, ref, minMatchL, currentMaxErrors, kValues, 1) // Stride 1
		for _, m := range chunkMatches {
			matches = append(matches, common.AnchorMatch{
				QueryStart:  m.QueryStart + chunkStart, // Adjust to queryRegion coordinates
				QueryEnd:    m.QueryEnd + chunkStart,
				RefStart:    m.RefStart,
				RefEnd:      m.RefEnd,
				Score:       m.Score,
				Identity:    m.Identity,
				Orientation: m.Orientation,
			})
		}
	}

	// Enhanced filtering (Python logic)
	if len(matches) == 0 {
		return []common.AnchorMatch{}
	}

	sort.SliceStable(matches, func(i, j int) bool { return matches[i].Score > matches[j].Score })

	var filtered []common.AnchorMatch
	excludedIndices := make(map[int]bool)

	for i := 0; i < len(matches); i++ {
		if excludedIndices[i] {
			continue
		}
		matchI := matches[i]
		filtered = append(filtered, matchI)

		qStartI, qEndI := matchI.QueryStart, matchI.QueryEnd
		lenI := (qEndI - qStartI) + 1

		overlapThresh := config.DefaultOverlapThreshold
		if matchI.Identity > 0.9 && lenI > 100 {
			overlapThresh = config.HighQualityOverlapThreshold
		}

		for j := 0; j < len(matches); j++ { // Python iterated all j, not just j > i
			if i == j || excludedIndices[j] {
				continue
			}
			matchJ := matches[j]
			qStartJ, qEndJ := matchJ.QueryStart, matchJ.QueryEnd

			overlapQStart := int(math.Max(float64(qStartI), float64(qStartJ)))
			overlapQEnd := int(math.Min(float64(qEndI), float64(qEndJ)))
			overlapLen := 0
			if overlapQEnd >= overlapQStart {
				overlapLen = overlapQEnd - overlapQStart + 1
			}

			minLen := int(math.Min(float64(lenI), float64((qEndJ-qStartJ)+1)))
			overlapRatio := 0.0
			if minLen > 0 {
				overlapRatio = float64(overlapLen) / float64(minLen)
			}

			if overlapRatio > overlapThresh {
				excludedIndices[j] = true
			}
		}
	}
	return filtered
}
