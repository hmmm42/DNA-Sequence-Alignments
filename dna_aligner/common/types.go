package common

// Segment represents a matched region between query and reference.
// QueryStart, QueryEnd, RefStart, RefEnd are 0-based inclusive.
type Segment struct {
	QueryStart int
	QueryEnd   int
	RefStart   int
	RefEnd     int
}

// AnchorMatch stores information about an extended k-mer match.
// Score and Identity are also stored.
// Orientation: 'f' for forward, 'r' for reverse (used internally).
type AnchorMatch struct {
	QueryStart  int // Inclusive
	QueryEnd    int // Inclusive
	RefStart    int // Inclusive
	RefEnd      int // Inclusive
	Score       float64
	Identity    float64
	Orientation rune // 'f' for forward, 'r' for reverse
}

// KmerMatch stores a simple k-mer exact match position.
type KmerMatch struct {
	QueryPos int
	RefPos   int
	Length   int
}

// Edge for the segment graph
type Edge struct {
	To     int     // Index of the target anchor/node
	Weight float64 // Score associated with transitioning to/including the target node
}
