package matching

import (
	"DNA-Sequence-Alignments/dna_aligner/common"
	"DNA-Sequence-Alignments/dna_aligner/config"
)

// FindExactMatches finds exact matches of length k between query and reference.
func FindExactMatches(query, ref string, k int) []common.KmerMatch {
	if k == 0 {
		k = config.DefaultK
	}
	if k <= 0 || k > len(ref) || k > len(query) { // Basic validation
		return []common.KmerMatch{}
	}

	refKmers := make(map[string][]int)
	for i := 0; i <= len(ref)-k; i++ {
		kmer := ref[i : i+k]
		refKmers[kmer] = append(refKmers[kmer], i)
	}

	var matches []common.KmerMatch
	for i := 0; i <= len(query)-k; i++ {
		kmer := query[i : i+k]
		if rPositions, found := refKmers[kmer]; found {
			for _, rPos := range rPositions {
				matches = append(matches, common.KmerMatch{QueryPos: i, RefPos: rPos, Length: k})
			}
		}
	}
	return matches
}
