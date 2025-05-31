package main

import "maps"

type State struct {
	len  int
	link int
	next map[byte]int
}

func NewState(len, link int) *State {
	return &State{len: len, link: link, next: make(map[byte]int)}
}

type SAM struct {
	last   int
	size   int
	states []*State
}

func NewSAM() *SAM {
	return &SAM{
		last:   0,
		size:   1,
		states: []*State{NewState(0, -1)},
	}
}

func BuildSAMByString(s string) *SAM {
	sam := NewSAM()
	for _, c := range s {
		sam.Extend(byte(c))
	}
	return sam
}

func (s *SAM) Extend(c byte) {
	p, cur := s.last, s.size
	s.size++
	s.states = append(s.states, NewState(s.states[p].len+1, -1))
	s.states[cur].len = s.states[p].len + 1

	for ; p != -1 && s.states[p].next[c] == 0; p = s.states[p].link {
		s.states[p].next[c] = cur
	}

	if p == -1 {
		s.states[cur].link = 0
	} else {
		q := s.states[p].next[c]
		if s.states[p].len+1 == s.states[q].len {
			s.states[cur].link = q
		} else {
			clone := s.size
			s.size++
			s.states = append(s.states, NewState(s.states[p].len+1, s.states[q].link))

			maps.Copy(s.states[clone].next, s.states[q].next)

			for ; p != -1 && s.states[p].next[c] == q; p = s.states[p].link {
				s.states[p].next[c] = clone
			}
			s.states[q].link = clone
			s.states[cur].link = clone
		}
	}
	s.last = cur
}

func (s *SAM) FindMaxMatch(query string, start int) int {
	maxLen := 0
	cur := 0
	for i := start; i < len(query); i++ {
		c := query[i]
		if _, ok := s.states[cur].next[c]; !ok {
			break
		}
		cur = s.states[cur].next[c]
		maxLen++
	}
	return maxLen

}
