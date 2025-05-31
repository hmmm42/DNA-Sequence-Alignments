package config

// K-mer matching parameters
const (
	DefaultK         = 10
	MinMatchLength   = 28
	DefaultMaxErrors = 5
	DefaultStride    = 2
)

// Extension parameters
const (
	ExtendMaxErrors      = 6
	MinIdentityThreshold = 0.74 // float64
)

// Anchor filtering parameters
const (
	DefaultOverlapThreshold     = 0.72 // float64
	HighQualityOverlapThreshold = 0.48 // float64
)

// Large region processing parameters
const (
	MaxSegmentSize       = 475
	LargeRegionChunkSize = 575
	LargeRegionOverlap   = 210
)
const StandardChunkOverlapRatio = 2.8 // float64

// Small region processing parameters
const (
	VerySmallRegionThreshold = 4
	SmallKForShortSegments   = 5
	MediumKForShortSegments  = 6
	LargeKForShortSegments   = 7
)

// Coverage parameters
const (
	SmallSegmentLength   = 275
	SamplePositionsCount = 25
)

// Merging parameters
const (
	AdjacentMergeMaxGap = 32
	FinalMergeMaxGap    = 22
)
const MaxGapRatioDifference = 0.55 // float64

// Adaptive parameters based on sequence properties
const (
	LowGCThreshold  = 0.40 // float64
	HighGCThreshold = 0.50 // float64
)

// Length-based parameters
const (
	ShortSeqThreshold = 3250
)

// K-mer size options based on sequence properties
var LowGCKValues = []int{8, 9, 10}
var MedGCKValues = []int{7, 8, 9}
var HighGCKValues = []int{6, 7, 8}

// Error tolerance based on sequence GC content
const (
	LowGCMaxErrors  = 4
	HighGCMaxErrors = 6
)

// K-mer options for different segment lengths
var VeryShortSegmentKValues = []int{4, 5}
var ShortSegmentKValues = []int{5, 6, 7}
var LongerSegmentKValues = []int{7, 8, 9}

// Additional error tolerance for boundary regions
const (
	BoundaryExtraErrors = 2
)

// Dataset specific optimizations (for very short sequences)
const (
	VeryShortSeqThreshold      = 2500
	VeryShortSeqMinMatchLength = 20
	VeryShortSeqStride         = 1
)

var VeryShortSeqKValues = []int{5, 6, 7}
