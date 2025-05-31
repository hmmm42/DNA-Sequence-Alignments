package main

import (
	"DNA-Sequence-Alignments/dna_aligner/aligner"
	"DNA-Sequence-Alignments/dna_aligner/common"
	"DNA-Sequence-Alignments/dna_aligner/io"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"
)

func formatSegmentsOutput(segments []common.Segment) string {
	var parts []string
	for _, seg := range segments {
		// Format: (query_start, query_end, ref_start, ref_end) - all inclusive
		parts = append(parts, fmt.Sprintf("(%d, %d, %d, %d)", seg.QueryStart, seg.QueryEnd+1, seg.RefStart, seg.RefEnd+1))
	}
	// Matches Python's str() output for a list of tuples: [(v1,v2,v3,v4), (v5,v6,v7,v8)]
	return "[" + strings.Join(parts, ", ") + "]"
}

func main() {
	// Determine data directory. Assumes 'data' is a subdirectory where the app is run,
	// or relative to the executable path.
	dataDir := "data" // Default to CWD/data
	exePath, err := os.Executable()
	if err == nil {
		dataDir = filepath.Join(filepath.Dir(exePath), "data")
	}
	// Fallback if executable path is weird or 'data' isn't there, try CWD/data again
	if _, statErr := os.Stat(dataDir); os.IsNotExist(statErr) {
		cwd, _ := os.Getwd()
		dataDir = filepath.Join(cwd, "data")
	}

	fmt.Printf("Using data directory: %s\n", dataDir)
	if _, statErr := os.Stat(dataDir); os.IsNotExist(statErr) {
		fmt.Printf("Error: Data directory '%s' not found. Please create it and place query/ref files.\n", dataDir)
		fmt.Println("Expected files like 'query1.txt', 'ref1.txt', etc., inside the data directory.")
		return
	}

	for i := 1; i <= 2; i++ { // Process dataset 1 and 2 as in Python's main
		queryFileName := fmt.Sprintf("query%d.txt", i)
		refFileName := fmt.Sprintf("ref%d.txt", i)
		queryFile := filepath.Join(dataDir, queryFileName)
		refFile := filepath.Join(dataDir, refFileName)
		outputFile := fmt.Sprintf("result%d.txt", i) // Output in CWD

		fmt.Printf("\nProcessing dataset %d (Query: %s, Ref: %s)...\n", i, queryFileName, refFileName)

		querySeq, err := io.ReadSequence(queryFile)
		if err != nil {
			fmt.Printf("Error reading query file '%s': %v\n", queryFile, err)
			continue
		}
		refSeq, err := io.ReadSequence(refFile)
		if err != nil {
			fmt.Printf("Error reading reference file '%s': %v\n", refFile, err)
			continue
		}

		fmt.Printf("Query length: %d, Reference length: %d\n", len(querySeq), len(refSeq))

		startTime := time.Now()

		// The aligner.FindAlignment function uses 0 for minMatchLenUser to trigger default/adaptive logic.
		alignmentResultSegments := aligner.FindAlignment(querySeq, refSeq, 0)

		duration := time.Since(startTime)
		fmt.Printf("Time taken for dataset %d: %.2f seconds\n", i, duration.Seconds())
		fmt.Printf("Found %d matching regions for dataset %d\n", len(alignmentResultSegments), i)

		outputString := formatSegmentsOutput(alignmentResultSegments)
		err = os.WriteFile(outputFile, []byte(outputString), 0644)
		if err != nil {
			fmt.Printf("Error writing output file '%s': %v\n", outputFile, err)
			continue
		}
		fmt.Printf("Results for dataset %d written to '%s'\n", i, outputFile)
	}
}
