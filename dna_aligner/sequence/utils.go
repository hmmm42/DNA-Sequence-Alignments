package sequence

import "strings"

// ReverseComplement returns the reverse complement of a DNA sequence.
func ReverseComplement(seq string) string {
	complement := map[rune]rune{
		'A': 'T', 'T': 'A',
		'C': 'G', 'G': 'C',
		'N': 'N', // Handle N if present
	}
	runes := []rune(seq)
	n := len(runes)
	var sb strings.Builder
	sb.Grow(n)
	for i := n - 1; i >= 0; i-- {
		base := runes[i]
		if compBase, ok := complement[base]; ok {
			sb.WriteRune(compBase)
		} else {
			sb.WriteRune('N') // Default for unknown bases
		}
	}
	return sb.String()
}

// CalculateGCContent calculates the GC content of a DNA sequence.
func CalculateGCContent(seq string) float64 {
	if len(seq) == 0 {
		return 0.0
	}
	gcCount := 0
	for _, base := range seq {
		if base == 'G' || base == 'C' || base == 'g' || base == 'c' { // Case-insensitive
			gcCount++
		}
	}
	return float64(gcCount) / float64(len(seq))
}
