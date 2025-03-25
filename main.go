package main

import (
	"os"
	"slices"
	"strings"
)

func init() {
	baseMapping = map[byte]byte{
		'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
	}
}

func readFile(f string) string {
	res, err := os.ReadFile(f)
	if err != nil {
		panic(err)
	}
	return strings.TrimSpace(string(res))
}

func reverseComplement(s string) string {
	bytes := []byte(s)
	slices.Reverse(bytes)
	for i, b := range bytes {
		bytes[i] = baseMapping[b]
	}
	return string(bytes)
}

func main() {
	ref := readFile("data/ref.txt")
	query := readFile("data/query.txt")

	res := analyzeDuplicates(query, ref)
	printDupIdentificationResults(res)
}
