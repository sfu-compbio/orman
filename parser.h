/* 786 
 *
 * Copyright (c) 2012, 2013, Simon Fraser University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the name of the Simon Fraser University nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author         : Ibrahim Numanagic
 * Email          : inumanag AT sfu DOT ca
 * Last Update    : 30. ix 2013.
 */

#ifndef PARSER_H__
#define PARSER_H__

#include "common.h"
#include "annotation.h"
#include "interval.h"
#include "partial.h"
using namespace std;

enum indel_type { INSERT, DELETE };
struct indel {
	int first, second;
	indel_type type;
	string insert;

	indel (void) {}
	indel (int f, int s, indel_type t) :
		first(f), second(s), type(t) {}
};

struct read_entry_key {
	int32_t  chr1, chr2;
	uint32_t pos1, pos2;

	int32_t  tlen;
	char     is_secondary;

	read_entry_key() {}
	read_entry_key(int32_t c1, uint32_t p1, int32_t c2, uint32_t p2, int32_t t, char i):
		chr1(c1), pos1(p1), chr2(c2), pos2(p2), tlen(t), is_secondary(i) {}
	bool operator< (const read_entry_key& x) const {
		if (chr1 != x.chr1) return chr1 < x.chr1;
		if (chr2 != x.chr2) return chr2 < x.chr2;
		if (pos1 != x.pos1) return pos1 < x.pos1;
		if (pos2 != x.pos2) return pos2 < x.pos2;
		if (tlen != x.tlen) return tlen < x.tlen;
		return is_secondary < x.is_secondary;
	}
};

struct read_entry_value {
	int64_t line1, line2;
	vector<interval> part1, part2;
	vector<indel> indel1, indel2;
};

struct read_pt_entry {
	PTsp partial;
	int chr1, chr2;
	int64_t line1, line2;
	int32_t start1, start2;
	vector<interval> part1, part2; 

	read_pt_entry() {}
	read_pt_entry(PTsp p, int c1, int64_t l1, int32_t s1, const vector<interval> &i1,
				          int c2, int64_t l2, int32_t s2, const vector<interval> &i2):
		partial(p), 
		chr1(c1), line1(l1), start1(s1), part1(i1),
		chr2(c2), line2(l2), start2(s2), part2(i2) {}
	bool operator< (const read_pt_entry& x) const {
		return (line1 < x.line1) || (line1 == x.line1 && line2 < x.line2);
	}
};

//extern set<PTs> PTs_set;
extern genome_annotation ga;
extern vector<struct read> reads;
extern vector<pair<int64_t, genome_annotation::gene*> > single_maps;
extern vector<int64_t> crappy_reads;
extern int read_length;

void parse_sam (const char *sam_file);

#endif // PARSER_H__