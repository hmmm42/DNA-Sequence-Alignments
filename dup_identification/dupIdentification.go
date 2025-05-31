package main

import (
	"fmt"
	"strings"
)

var baseMapping map[byte]byte

type Duplicate struct {
	QueryStart int
	RefStart   int
	Length     int
	Count      int
	IsInverted bool
}

func analyzeDuplicates(query, ref string) []Duplicate {

	// 生成反向互补序列
	invRef := reverseComplement(ref)

	// 构建正向和反向引用序列的后缀自动机
	refSAM := BuildSAMByString(ref)
	invRefSAM := BuildSAMByString(invRef)

	// 存储结果
	var duplicates []Duplicate

	// 预计算查询序列每个位置的最优匹配
	type matchInfo struct {
		maxLength int
		isInv     bool
	}

	// 为查询序列的每个位置计算最大匹配长度和匹配类型
	matchData := make([]matchInfo, len(query))
	for pos := 0; pos < len(query); pos++ {
		// 计算正向和反向的最大匹配长度
		forwardLen := refSAM.FindMaxMatch(query, pos)
		reverseLen := invRefSAM.FindMaxMatch(query, pos)

		// 决定使用哪种匹配
		isInverted := reverseLen > forwardLen || (reverseLen == forwardLen && reverseLen > 0)
		bestLen := forwardLen
		if isInverted {
			bestLen = reverseLen
		}

		matchData[pos] = matchInfo{
			maxLength: bestLen,
			isInv:     isInverted,
		}
	}

	// 处理查询序列，寻找所有可能的重复
	position := 0
	for position < len(query) {
		currentMatch := matchData[position]

		// 如果当前位置没有匹配，移动到下一个位置
		if currentMatch.maxLength == 0 {
			position++
			continue
		}

		// 获取重复单元的长度和类型
		unitLength := currentMatch.maxLength
		unitInverted := currentMatch.isInv

		// 提取重复单元
		repeatUnit := query[position : position+unitLength]

		// 计算连续重复次数
		repeatCount := 1
		nextStart := position + unitLength

		for nextStart+unitLength <= len(query) {
			// 检查下一个单元是否符合重复条件
			if query[nextStart:nextStart+unitLength] != repeatUnit ||
				matchData[nextStart].maxLength < unitLength ||
				matchData[nextStart].isInv != unitInverted {
				break
			}
			repeatCount++
			nextStart += unitLength
		}

		// 找到重复单元在参考序列中的起始位置
		refUnit := repeatUnit
		if unitInverted {
			refUnit = reverseComplement(repeatUnit)
		}

		// 在参考序列中查找第一次出现的位置
		refPosition := strings.Index(ref, refUnit)

		// 记录结果
		duplicates = append(duplicates, Duplicate{
			QueryStart: position,
			RefStart:   refPosition,
			Length:     unitLength,
			Count:      repeatCount,
			IsInverted: unitInverted,
		})

		// 跳到下一个未处理的位置
		position = nextStart
	}

	return duplicates
}

func printDupIdentificationResults(results []Duplicate) {
	fmt.Println("Duplicate Identification Results")
	fmt.Println("|   Pos in Ref   |   Repeat Size   |   Repeat Count   |   Inverse   |")
	fmt.Println("|----------------|-----------------|------------------|-------------|")

	for _, entry := range results {
		invertedStr := "Yes"
		if !entry.IsInverted {
			invertedStr = "No"
		}
		fmt.Printf("|   %-12d |   %-13d |   %-14d |   %-9s |\n",
			entry.RefStart, entry.Length, entry.Count, invertedStr)
	}
}
