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

#ifndef PARTIAL_H__
#define PARTIAL_H__

#include "common.h"
#include "annotation.h"

using namespace std;

struct PT_single {
	// TODO gene only for now
	genome_annotation::transcript *transcript;
	string signature;
	int length;
	int weight;
	uint32_t start;

	PT_single (void) {}
	PT_single (genome_annotation::transcript *t, const string &s, int l, int w, uint32_t st) : 
		transcript(t), signature(s), length(l), weight(w), start(st) {}

	bool operator== (const PT_single &s) const {
		return (transcript == s.transcript && signature == s.signature);
	}
	bool operator< (const PT_single &s) const {
		return (transcript < s.transcript || (transcript == s.transcript && signature < s.signature));
	}

	genome_annotation::gene *get_gene(void) const {
		return transcript->gene;
	}
};
typedef PT_single 			PTs;
typedef const PT_single* 	PTsp;
typedef PTsp PT;

 uint32_t get_absolute_position (const PT &p, int k);
 uint32_t get_gene_position (const PT &p, int k);
 uint32_t g2G ( genome_annotation::transcript *t, uint32_t k) ;
 int get_single_coverage(const PT &p, int k);
 int get_single_coverage(uint32_t j);
 char is_multimap (uint32_t j);

struct read {
	struct read_entry {
		PT partial;
		pair<uint32_t, uint32_t> partial_start;
		pair<int, int> 			 line;
		
		read_entry (void) {}
		read_entry (const PT &p, const pair<int, int> &l, const pair<uint32_t, uint32_t> &s) : 
			partial(p), line(l), partial_start(s) {}
	};
	vector<read_entry> entries;
};

#endif // PARTIAL_H__

