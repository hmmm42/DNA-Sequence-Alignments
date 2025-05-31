package regions

import (
	"DNA-Sequence-Alignments/dna_aligner/common"
	"DNA-Sequence-Alignments/dna_aligner/config"
	"fmt"
	"math"
	"math/rand"
	"sort"
	"time"
)

func init() {
	rand.Seed(time.Now().UnixNano())
}

// ResolveOverlaps resolves overlapping segments in query coordinates by keeping the longer segment.
// Input segments are sorted by QueryStart. If not, this function will sort them.
func ResolveOverlaps(segments []common.Segment) []common.Segment {
	if len(segments) == 0 {
		return []common.Segment{}
	}

	sortedSegs := make([]common.Segment, len(segments))
	copy(sortedSegs, segments)
	sort.Slice(sortedSegs, func(i, j int) bool {
		if sortedSegs[i].QueryStart != sortedSegs[j].QueryStart {
			return sortedSegs[i].QueryStart < sortedSegs[j].QueryStart
		}
		// If query starts are same, prioritize longer or one with smaller ref start
		lenI := sortedSegs[i].QueryEnd - sortedSegs[i].QueryStart
		lenJ := sortedSegs[j].QueryEnd - sortedSegs[j].QueryStart
		if lenI != lenJ {
			return lenI > lenJ
		} // Longer first
		return sortedSegs[i].RefStart < sortedSegs[j].RefStart
	})

	var result []common.Segment
	if len(sortedSegs) > 0 {
		result = append(result, sortedSegs[0])
		for i := 1; i < len(sortedSegs); i++ {
			curr := sortedSegs[i]
			prev := &result[len(result)-1] // Get pointer to last segment in result

			if curr.QueryStart <= prev.QueryEnd { // Overlap
				currLen := (curr.QueryEnd - curr.QueryStart) + 1
				prevLen := (prev.QueryEnd - prev.QueryStart) + 1
				if currLen > prevLen { // Keep longer
					*prev = curr
				} // Else, prev (which is longer or equal and came first) is kept
			} else { // No overlap
				result = append(result, curr)
			}
		}
	}
	return result
}

// isSegmentOriginal checks if a segment was part of the initial set.
func isSegmentOriginal(segment common.Segment, originalSegments []common.Segment) bool {
	for _, origSeg := range originalSegments {
		if segment.QueryStart == origSeg.QueryStart &&
			segment.QueryEnd == origSeg.QueryEnd &&
			segment.RefStart == origSeg.RefStart &&
			segment.RefEnd == origSeg.RefEnd {
			return true
		}
	}
	return false
}

// EnsureCompleteCoverage ensures the entire query is covered by finding matches for uncovered regions.
func EnsureCompleteCoverage(query, ref string, initialSegments []common.Segment) []common.Segment {
	queryLen := len(query)
	refLen := len(ref)
	if queryLen == 0 {
		return []common.Segment{}
	}

	// initialSegments should be sorted for FindUncoveredRegions.
	// Let's make a mutable copy and sort it.
	currentCoverageSegments := make([]common.Segment, len(initialSegments))
	copy(currentCoverageSegments, initialSegments)
	sort.Slice(currentCoverageSegments, func(i, j int) bool {
		return currentCoverageSegments[i].QueryStart < currentCoverageSegments[j].QueryStart
	})

	uncovered := FindUncoveredRegions(queryLen, currentCoverageSegments)

	if len(uncovered) == 0 {
		fmt.Println("Query already has complete coverage based on initial segments.")
		return currentCoverageSegments // Already sorted and presumably non-overlapping if initialSegments were clean
	}
	fmt.Printf("Found %d uncovered regions in query\n", len(uncovered))

	newlyFoundSegments := []common.Segment{}

	for i, regionCoords := range uncovered {
		qStart, qEnd := regionCoords[0], regionCoords[1]
		regionActualLen := (qEnd - qStart) + 1
		fmt.Printf("Processing uncovered region %d/%d: query pos %d-%d (length: %d)\n",
			i+1, len(uncovered), qStart, qEnd, regionActualLen)

		if regionActualLen < config.VerySmallRegionThreshold {
			fmt.Printf("  Skipping very small region (length: %d)\n", regionActualLen)
			continue
		}

		queryRegionStr := query[qStart : qEnd+1]
		var regionMatches []common.AnchorMatch // Relative coordinates

		if regionActualLen > 1000 { // Python's threshold for "large region" specific handling
			fmt.Println("  Large region detected, using divide-and-conquer approach")
			regionMatches = FindMatchesInLargeRegion(queryRegionStr, ref, 500, 0, 0) // 500 is from python's call
		} else {
			regionMatches = FindMatchesInRegion(queryRegionStr, ref, 0, 0)
		}

		if len(regionMatches) > 0 {
			fmt.Printf("  Found %d potential matches for this region\n", len(regionMatches))
			sort.Slice(regionMatches, func(i, j int) bool { // Sort by score
				return regionMatches[i].Score > regionMatches[j].Score
			})

			// Add non-overlapping (with current new finds for this region) matches
			tempAddedForThisRegion := []common.Segment{}
			for _, match := range regionMatches {
				absQStart := match.QueryStart + qStart
				absQEnd := match.QueryEnd + qStart

				// Check overlap with segments already added *for this specific uncovered region processing iteration*
				isOverlappingTemp := false
				for _, tempSeg := range tempAddedForThisRegion {
					if math.Max(float64(tempSeg.QueryStart), float64(absQStart)) <= math.Min(float64(tempSeg.QueryEnd), float64(absQEnd)) {
						isOverlappingTemp = true
						break
					}
				}
				if !isOverlappingTemp {
					segToAdd := common.Segment{
						QueryStart: absQStart, QueryEnd: absQEnd,
						RefStart: match.RefStart, RefEnd: match.RefEnd,
					}
					tempAddedForThisRegion = append(tempAddedForThisRegion, segToAdd)
					newlyFoundSegments = append(newlyFoundSegments, segToAdd)
				}
			}
		} else { // No matches found for region, Python's fallback logic
			fmt.Printf("  No matches found for region. Creating fallback segments.\n")
			if regionActualLen > config.MaxSegmentSize {
				chunkStep := config.SmallSegmentLength
				for chunkOffset := 0; chunkOffset < regionActualLen; chunkOffset += chunkStep {
					curChunkStartInRegion := chunkOffset
					curChunkEndInRegion := int(math.Min(float64(chunkOffset+chunkStep), float64(regionActualLen)))
					if curChunkEndInRegion <= curChunkStartInRegion {
						continue
					}

					absQChunkStart := qStart + curChunkStartInRegion
					absQChunkEnd := qStart + curChunkEndInRegion - 1
					chunkStr := query[absQChunkStart : absQChunkEnd+1]
					chunkActualLen := len(chunkStr)

					bestRStart, bestScore := 0, -1.0
					sampleStep := int(math.Max(1, float64(refLen/config.SamplePositionsCount)))
					if sampleStep == 0 {
						sampleStep = 1
					}

					for rPos := 0; rPos <= refLen-chunkActualLen; rPos += sampleStep {
						refChunkStr := ref[rPos : rPos+chunkActualLen]
						matchesCount := 0
						for k := 0; k < chunkActualLen; k++ {
							if chunkStr[k] == refChunkStr[k] {
								matchesCount++
							}
						}
						score := float64(matchesCount) / float64(chunkActualLen)
						if score > bestScore {
							bestScore, bestRStart = score, rPos
						}
					}
					newlyFoundSegments = append(newlyFoundSegments, common.Segment{
						QueryStart: absQChunkStart, QueryEnd: absQChunkEnd,
						RefStart: bestRStart, RefEnd: bestRStart + chunkActualLen - 1,
					})
				}
			} else { // Smaller region, no matches, add one fallback segment
				rMapStart := qStart % refLen
				rMapEnd := rMapStart + regionActualLen - 1
				if rMapEnd >= refLen { // Ensure it fits
					rMapEnd = refLen - 1
					if rMapStart > rMapEnd && regionActualLen <= refLen {
						rMapStart = 0
					} // Adjust start if possible
					if regionActualLen > refLen {
						rMapStart = 0
					} // query region longer than ref
				}
				newlyFoundSegments = append(newlyFoundSegments, common.Segment{
					QueryStart: qStart, QueryEnd: qEnd, RefStart: rMapStart, RefEnd: rMapEnd,
				})
			}
		}
	}

	// Combine initial segments with newly found ones for uncovered regions
	allSegments := append(currentCoverageSegments, newlyFoundSegments...)
	sort.Slice(allSegments, func(i, j int) bool { // Sort all before resolving overlaps
		return allSegments[i].QueryStart < allSegments[j].QueryStart
	})

	// Python's complex non-overlapping resolution:
	// "Choose the better segment (prefer original segments, then longer ones)"
	var resolvedWithPreference []common.Segment
	if len(allSegments) > 0 {
		current := allSegments[0]
		for i := 1; i < len(allSegments); i++ {
			next := allSegments[i]
			if next.QueryStart <= current.QueryEnd { // Overlap
				isCurrentOrig := isSegmentOriginal(current, initialSegments)
				isNextOrig := isSegmentOriginal(next, initialSegments)

				if isNextOrig && !isCurrentOrig {
					current = next
				} else if !isNextOrig && isCurrentOrig {
					// keep current
				} else { // Both original or both new: keep longer
					if (next.QueryEnd - next.QueryStart) > (current.QueryEnd - current.QueryStart) {
						current = next
					}
				}
			} else {
				resolvedWithPreference = append(resolvedWithPreference, current)
				current = next
			}
		}
		resolvedWithPreference = append(resolvedWithPreference, current)
	}

	// Final check for uncovered regions and fill them if any (Python's final fallback)
	finalUncovered := FindUncoveredRegions(queryLen, resolvedWithPreference)
	if len(finalUncovered) > 0 {
		fmt.Printf("Warning: %d regions still uncovered. Adding final fallback segments.\n", len(finalUncovered))
		for _, reg := range finalUncovered {
			s, e := reg[0], reg[1]
			rLen := (e - s) + 1
			rS := s % refLen
			rE := rS + rLen - 1
			if rE >= refLen {
				rE = refLen - 1
				if rS > rE && rLen <= refLen {
					rS = 0
				} else if rLen > refLen {
					rS = 0
				}
			}
			resolvedWithPreference = append(resolvedWithPreference, common.Segment{QueryStart: s, QueryEnd: e, RefStart: rS, RefEnd: rE})
		}
		// Re-sort and resolve all overlaps finally
		return ResolveOverlaps(resolvedWithPreference)
	}

	return resolvedWithPreference
}
