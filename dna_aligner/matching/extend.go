package matching

import (
	"DNA-Sequence-Alignments/dna_aligner/common"
	"DNA-Sequence-Alignments/dna_aligner/config"
	"math"
)

// ExtendMatch extends a k-mer match to a longer anchor with error tolerance.
// Returns an AnchorMatch with inclusive coordinates if a valid extension is found, otherwise nil.
func ExtendMatch(query, ref string, qStartKmer, rStartKmer, k int, minMatchLengthUser int, maxErrorsUser int) *common.AnchorMatch {
	minMatchLen := minMatchLengthUser
	if minMatchLen == 0 {
		minMatchLen = config.MinMatchLength
	}
	maxErrors := maxErrorsUser
	if maxErrors == 0 {
		maxErrors = config.ExtendMaxErrors
	}

	// Initial state from k-mer (exclusive ends for loop variables)
	qCurrentFwd, rCurrentFwd := qStartKmer+k, rStartKmer+k
	totalMatches := k
	errorsFwd := 0

	// Extend forward
	for qCurrentFwd < len(query) && rCurrentFwd < len(ref) && errorsFwd <= maxErrors {
		if query[qCurrentFwd] == ref[rCurrentFwd] {
			qCurrentFwd++
			rCurrentFwd++
			totalMatches++
		} else {
			indelMatchFound := false
			// Try insertion in query (gap in ref)
			for ins := 1; ins <= 2; ins++ {
				if qCurrentFwd+ins < len(query) && rCurrentFwd < len(ref) && query[qCurrentFwd+ins] == ref[rCurrentFwd] {
					qCurrentFwd += (ins + 1)
					rCurrentFwd++
					errorsFwd++    // Indel costs 1 error
					totalMatches++ // The matching base after indel
					indelMatchFound = true
					break
				}
			}
			if indelMatchFound {
				continue
			}

			// Try insertion in reference (gap in query)
			for ins := 1; ins <= 2; ins++ {
				if rCurrentFwd+ins < len(ref) && qCurrentFwd < len(query) && query[qCurrentFwd] == ref[rCurrentFwd+ins] {
					qCurrentFwd++
					rCurrentFwd += (ins + 1)
					errorsFwd++
					totalMatches++
					indelMatchFound = true
					break
				}
			}
			if indelMatchFound {
				continue
			}

			// Mismatch
			qCurrentFwd++
			rCurrentFwd++
			errorsFwd++
		}
	}
	// qCurrentFwd, rCurrentFwd are now exclusive ends for the forward extended part.

	// Extend backward (inclusive starts for loop variables)
	qCurrentBwd, rCurrentBwd := qStartKmer-1, rStartKmer-1
	errorsBwd := 0 // Reset error count for backward extension, as in Python

	for qCurrentBwd >= 0 && rCurrentBwd >= 0 && errorsBwd <= maxErrors {
		if query[qCurrentBwd] == ref[rCurrentBwd] {
			qCurrentBwd--
			rCurrentBwd--
			totalMatches++
		} else {
			indelMatchFound := false
			// Try insertion in query (gap in ref) - looking backward
			for ins := 1; ins <= 2; ins++ {
				if qCurrentBwd-ins >= 0 && rCurrentBwd >= 0 && query[qCurrentBwd-ins] == ref[rCurrentBwd] {
					qCurrentBwd -= (ins + 1)
					rCurrentBwd--
					errorsBwd++
					totalMatches++
					indelMatchFound = true
					break
				}
			}
			if indelMatchFound {
				continue
			}

			// Try insertion in reference (gap in query) - looking backward
			for ins := 1; ins <= 2; ins++ {
				if rCurrentBwd-ins >= 0 && qCurrentBwd >= 0 && query[qCurrentBwd] == ref[rCurrentBwd-ins] {
					qCurrentBwd--
					rCurrentBwd -= (ins + 1)
					errorsBwd++
					totalMatches++
					indelMatchFound = true
					break
				}
			}
			if indelMatchFound {
				continue
			}

			// Mismatch
			qCurrentBwd--
			rCurrentBwd--
			errorsBwd++
		}
	}
	// Final inclusive start positions
	finalQStart := qCurrentBwd + 1
	finalRStart := rCurrentBwd + 1

	// Final exclusive end positions are qCurrentFwd, rCurrentFwd
	finalQEndExclusive := qCurrentFwd
	finalREndExclusive := rCurrentFwd

	matchLength := finalQEndExclusive - finalQStart
	if matchLength <= 0 {
		return nil
	}

	identity := float64(totalMatches) / float64(matchLength)

	// Context adjustment from Python (applied to identity)
	if matchLength > 50 {
		contextSize := int(math.Min(20, float64(matchLength/4)))
		if finalQStart >= contextSize && finalRStart >= contextSize {
			leftQ := query[finalQStart-contextSize : finalQStart]
			leftR := ref[finalRStart-contextSize : finalRStart]
			leftMatchesCount := 0
			for i := 0; i < contextSize; i++ {
				if leftQ[i] == leftR[i] {
					leftMatchesCount++
				}
			}
			identity = (identity*float64(matchLength) + float64(leftMatchesCount)*0.5) / (float64(matchLength) + float64(contextSize)*0.5)
		}
	}

	// Python's scoring used the `errors` variable which, due to scope, was the one from backward extension (`errorsBwd`).
	scoreRelevantErrors := errorsBwd

	if matchLength >= minMatchLen && identity >= config.MinIdentityThreshold {
		score := float64(matchLength) * identity * (1.0 - 0.05*float64(scoreRelevantErrors))
		return &common.AnchorMatch{
			QueryStart: finalQStart,
			QueryEnd:   finalQEndExclusive - 1, // Store inclusive end
			RefStart:   finalRStart,
			RefEnd:     finalREndExclusive - 1, // Store inclusive end
			Score:      score,
			Identity:   identity,
		}
	}
	return nil
}
