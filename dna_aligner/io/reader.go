package io

import (
	"os"
	"strings"
)

// ReadSequence reads a DNA sequence from a file.
func ReadSequence(filePath string) (string, error) {
	data, err := os.ReadFile(filePath)
	if err != nil {
		return "", err
	}
	return strings.TrimSpace(string(data)), nil
}
