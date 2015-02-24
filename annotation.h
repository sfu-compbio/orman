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

#ifndef ANNOTATION_H__
#define ANNOTATION_H__

#include "common.h"
#include "interval.h"

class genome_annotation {
public:
	struct transcript;
	struct gene;

	struct exon {
		uint32_t start,
				 end;
		int id;
		string sid;
		struct transcript *transcript;

		exon (void) {}
		exon (uint32_t s, uint32_t e, const string &i) :
			start(s), end(e), sid(i) {
		}
		bool operator< (const exon &e) const {
			return start < e.start;
		}
	};
	struct transcript {
		int id;
		string name;
		vector<exon> exons;
		struct gene *gene;
		uint32_t length_;

		transcript (void) {}

		uint32_t length (void) const {
			return length_;
		}

		int chromosome (void) const {
			return gene->chromosome;
		}

		int position (int eid, int k) const {
			// end exclusive
			assert(eid < exons.size());
			for (int i = eid; i < exons.size(); i++) {
				if (k <= exons[i].end - exons[i].start + 1)
					return exons[i].start + k;
				else 
					k -= exons[i].end - exons[i].start + 1;
			}
			throw k;
		}
	};
	struct gene {
		int chromosome;
		map<string, transcript> transcripts;
		string name;

		gene (void) {}
	};

private:
	map<string, gene> genes;
	map<string, int>  chromosomes;
	map<int, uint32_t> chromosome_offset;

	interval_tree<exon*> exon_tree;

	int transcript_count;

public:
	genome_annotation (void) {

	}

	void parse_gtf (const char *gtf_file);
	int get_chromosome (const string &s) {
		map<string, int>::iterator it = chromosomes.find(s);
		if (it != chromosomes.end()) return it->second;
		else return -1;
	}

	void set_chromosome_offset(int chr, uint32_t x) {
		chromosome_offset[chr] = x;
	}
	uint32_t get_chromosome_offset(int chr) {
		return chromosome_offset[chr];
	}

	void get_exons (const interval &i, vector<genome_annotation::exon*> &exons);
	gene *find_gene (uint32_t pos, int chr);
	int get_transcript_count (void) const { return transcript_count; }

	void gt(map<string,transcript*> &m) {
		foreach(g,genes) {
			foreach(t,g->second.transcripts)
				m[t->first]=&t->second;
		}
	}
};


#endif // COMMON_H__


