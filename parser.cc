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

#include "parser.h"
#include <cstdarg>
using namespace std;

typedef genome_annotation::transcript transcript;
typedef genome_annotation::exon exon;
static string S (const char* f, ...) {
	char bf[MAX_BUFFER];
	va_list args;
	va_start(args, f);
	vsprintf(bf, f, args);
	va_end(args);
	return string(bf);
}


static unsigned short * covx;
static char *multiregion;
static uint32_t genome_len;
set<PTs> PTs_set;
vector<struct read> reads;
vector<pair<int64_t, genome_annotation::gene*> > single_maps;
vector<int64_t> crappy_reads;

uint32_t get_absolute_position (const PT &p, int k) {
	return 
		ga.get_chromosome_offset(p->transcript->chromosome()) +  // chr
		p->transcript->position(p->start, k); // position in exon ...
}

uint32_t get_gene_position (const PT &p, int k) {
	uint32_t e = k;
	for (int i = 0; i < p->start; i++)
		e += p->transcript->exons[i].end - p->transcript->exons[i].start + 1;
	assert(e < p->transcript->length());
	return e;
}

uint32_t g2G (genome_annotation::transcript *t, uint32_t k) {
 	assert (k < t->length());
	return ga.get_chromosome_offset(t->chromosome()) + 
		   t->position(0, k);
}

int get_single_coverage(const PT &p, int k) {
	return covx[get_absolute_position(p, k)];
}

int get_single_coverage(uint32_t j) {
	return covx[j];
}

char is_multimap (uint32_t j) {
	return multiregion[j];
}

void increase_coverage(int chr, int p1, int p2) {
	assert(p1 <= p2);
	uint32_t x = ga.get_chromosome_offset(chr);
	assert(x + p2 < genome_len + 1);
	for (int i = p1; i < p2; i++) 
		covx[x + i]++;	
}
/*****************************************************************************/

bool parse_cigar (uint32_t start_pos, const char *cigar, const char *read, vector<interval> &result, vector<indel> &indels) {
	int num = 0;

	if (strcmp(cigar, "*") == 0)
		return 0;

	while (*cigar) {
		if (isdigit(*cigar))
			num = 10 * num + (*cigar - '0');
		else {
			if (*cigar == 'M') {
				result.push_back(interval(start_pos, start_pos + num - 1)); // inclusive!!!
				start_pos += num;
				read += num;
			}
			else if (*cigar == 'D') { // delete FROM REFERENCE GENOME
				indels.push_back(indel(start_pos, start_pos + num, DELETE));
			//	start_pos += num;
			}
			else if (*cigar == 'N') // TODO check is the same as DEL
				start_pos += num;
			else if (*cigar == 'I') { // insert IN READ 
				indel x(start_pos, start_pos + num, INSERT);
				x.insert = string(read, num);
				indels.push_back(x);
				read += num;
			}
			else if (*cigar == 'S') {
				read += num;
			}
			else { // H O = X
				// E("CIGAR feature %c not implemented!", *cigar);
				return 0;
			}
			num = 0;
		}
		cigar++;
	}

	return 1;
}

#ifdef LOGIFY
static string get_partial_transcript_error;
static string _sline2, _sline1;
#endif

bool is_valid_read (const vector<interval> &parts, const vector<exon*> &exons, int &starting_position) {
	#ifdef LOGIFY
	get_partial_transcript_error = "";
	#endif	

	if (exons.size() == 0 || parts.size() == 0) 
		return 0;

	#ifdef LOGIFY
	exon *pt1 = exons[0];
	get_partial_transcript_error += S("%s:ERR=", string(pt1->transcript->gene->name /*+ "." + pt1->transcript->name*/).c_str());
	#endif

	// disregard reads with mapping chunk <4
	foreach (pi, parts)
		if (pi->second - pi->first + 1 < 4) {
			#ifdef LOGIFY
			get_partial_transcript_error += S("part_sz<4;");
			#endif
			return 0;
		}
	// disregard reads with introns and internal exon skips
	int pi = 0;
	foreach (ei, exons) {
		exon *e = *ei;
		// disregard introns
		// only allow +- 4 introns at the beginning
		if (parts[pi].first < e->start && !(e->start - parts[pi].first < 4)) { 
			#ifdef LOGIFY
			get_partial_transcript_error += S("intronS_exon_%d;", e->id);
			#endif
			return 0; 
		}
		// find internal skips
		int prev = parts[pi].first;
		while (pi < parts.size() && ((parts[pi].second <= e->end) || (parts[pi].second - e->end < 4))) {
			// internal exon skip
			if (prev < parts[pi].first) {
				#ifdef LOGIFY
				get_partial_transcript_error += S("inexon_skip_%d;", e->id);
				#endif
				return 0;
			}
			prev = parts[pi].second;
			pi++;
		}
		// again, no introns
		// allow +-4 at the end as well
		if (pi < parts.size () && parts[pi].first <= e->end && !(parts[pi].second - e->end < 4)) {
			#ifdef LOGIFY
			get_partial_transcript_error += S("intronE_exon_%d;", e->id);
			#endif
			return 0;
		}
	}
	starting_position = parts[0].first - exons[0]->start;
	if (starting_position < 0) 
		starting_position = 0; // +-4 mapper adjustment; larger cases are already discarded 
	return 1;
}

int get_single_weight (transcript *t, const vector<exon*> &exons, const vector<indel> &indels) {
	int weight = 1;
	int indel_idx = 0;
	for (int i = 0; i < exons.size(); i++) {
		exon *e = exons[i];
		// penalize missing exons by their length
		if (i) for (int j = exons[i - 1]->id + 1; j < exons[i]->id; j++)
			weight += 100 * (t->exons[j].end - t->exons[j].start + 1);
		// add exons
		// indels have to be in exons; otherwise we don't care
		// also they have to be +- 6 from sides
		while (indel_idx < indels.size() && indels[indel_idx].first < e->start + 6)
			indel_idx++;
		// we have indel!
		while (indel_idx < indels.size() && indels[indel_idx].second < e->end - 6) {
			// penalize indels
			weight += 10000;
			indel_idx++;
		}
	}
	return weight;
}

bool exonptr_sort (exon const *lhs, exon const *rhs) {
    return (*lhs) < (*rhs);
}
bool indel_sort (const indel &lhs, const indel &rhs) {
    return lhs.first < rhs.first;
}

PTsp get_partial_transcript (transcript *t, 
	const vector<interval> &part1, const vector<exon*> &exon1, const vector<indel> &indel1, int &start1,
	const vector<interval> &part2, const vector<exon*> &exon2, const vector<indel> &indel2, int &start2) 
{	
	bool v1 = is_valid_read(part1, exon1, start1);
	bool v2 = is_valid_read(part2, exon2, start2);
	if (!v1) return 0;
	int weight = get_single_weight(t, exon1, indel1) + 
				 v2 ? get_single_weight(t, exon2, indel2) : 1000000;
	if (!v2) start2 = 0;

	// indels, exons are sorted; detect them
	string signature = "";

	vector<exon*> exons;
	exons.insert(exons.end(), exon1.begin(), exon1.end());
	if (v2) exons.insert(exons.end(), exon2.begin(), exon2.end());
	sort(exons.begin(), exons.end(), exonptr_sort);

	vector<indel> indels;
	indels.insert(indels.end(), indel1.begin(), indel1.end());
	if (v2) indels.insert(indels.end(), indel2.begin(), indel2.end());
	sort(indels.begin(), indels.end(), indel_sort);

	int length = 0;
	int indel_idx = 0;
	for (int i = 0; i < exons.size(); i++) {
		exon *e = exons[i];
		if (i && e->id == exons[i - 1]->id) continue; // duplicates
		length += e->end - e->start + 1;
		// add exons
		signature += "e" + e->sid; //numtostr(e->id);
		// indels have to be in exons; otherwise we don't care
		// also they have to be +- 6 from sides
		while (indel_idx < indels.size() && indels[indel_idx].first < e->start + 6)
			indel_idx++;
		// we have indel!
		while (indel_idx < indels.size() && indels[indel_idx].second < e->end - 6) {
			signature += "[";
			if (indels[indel_idx].type == INSERT) { // reference insert
				signature += "I." + indels[indel_idx].insert;
			//	length += indels[indel_idx].second - indels[indel_idx].first + 1;
			}
			else { // reference delete
				signature += "D";
			//	length += indels[indel_idx].second - indels[indel_idx].first + 1;
			}
			signature += "," + numtostr(indels[indel_idx].first - e->start) // start
						  + "," + numtostr(indels[indel_idx].second - indels[indel_idx].first) // length
						  + "]";
			indel_idx++;
		}
	}

	PTs pt(t, signature, length, weight, exons[0]->id);
	pair<set<PTs>::iterator, bool> pti = PTs_set.insert(pt);
	return &(* pti.first);
}

void makecand (vector<interval> parts, int chromosome, map<transcript*, vector<exon*> > &candidates) {
	if (chromosome == -1) return;
	if (parts.size() == 0)
		return;
	if (parts[0].second - parts[0].first + 1 < 4)
		parts = vector<interval>(parts.begin() + 1, parts.end());
	if (parts[parts.size() - 1].second - parts[parts.size() - 1].first + 1 < 4)
		parts.pop_back();

	foreach (pi, parts) {
		vector<exon*> exons;
		ga.get_exons(*pi, exons);
		foreach (ei, exons) 
			if ((*ei)->transcript->gene->chromosome == chromosome) {
				vector<exon*> &vx = candidates[(*ei)->transcript];
				if (vx.size() == 0 || vx[vx.size() - 1] != *ei) // no repetitions!
					vx.push_back(*ei);
			}
	}
}

void parse_read (const read_entry_key &rk, const read_entry_value &rv, set<read_pt_entry> &result) {
	map<transcript*, vector<exon*> > candidates1, candidates2;
	makecand(rv.part1, rk.chr1, candidates1);
	makecand(rv.part2, rk.chr2, candidates2);
	int start1, start2;
	vector<exon*> dummy;
	foreach (ci1, candidates1) {
		bool no_null = false;
		// same transcripts only
		foreach (ci2, candidates2) if (ci1->first == ci2->first) {
			PTsp p = get_partial_transcript(ci1->first, 
				rv.part1, ci1->second, rv.indel1, start1,
				rv.part2, ci2->second, rv.indel2, start2
			);
			if (p) result.insert(read_pt_entry(p, 
						rk.chr1, rv.line1, start1, rv.part1,
						rk.chr2, rv.line2, start2, rv.part2));		
			if (p) no_null = true;
		}
		if (!no_null) {
			PTsp p = get_partial_transcript(ci1->first, 
				rv.part1, ci1->second, rv.indel1, start1,
				rv.part2, dummy, rv.indel2, start2
			);
			if (p) result.insert(read_pt_entry(p, 
						rk.chr1, rv.line1, start1, rv.part1,
						rk.chr2, rv.line2, start2, rv.part2));
		}
	}
}

void parse_sam (const char *sam_file) {
	FILE *fi = fopen(sam_file, "r");

	// obtain the file size
	fseek(fi, 0, SEEK_END);
	int64_t f_size = ftell(fi);
	fseek(fi, 0, SEEK_SET);

	// UGLY! but works for now ...
	multimap<read_entry_key, read_entry_value> idx;

	char *buffer = new char[MAX_BUFFER];
	char *sam_name = new char[MAX_BUFFER],
		 *sam_rname = new char[MAX_BUFFER],
		 *sam_cigar = new char[MAX_BUFFER],
		 *sam_read = new char[MAX_BUFFER],
		 *sam_rnext = new char[MAX_BUFFER];
	uint32_t sam_flag, sam_pos;
	int32_t sam_pnext, sam_tlen;
	uint8_t sam_mapq;

	int64_t line = 0;
	int read_id = 0;
	string prev_name = "";

	while (1) {
		fgets(buffer, MAX_BUFFER, fi);
		if (!feof(fi)) {
			if (buffer[0] == '@') {
				if (strlen(buffer) > 3 && buffer[1] == 'S' && buffer[2] == 'Q') {
					sscanf(buffer, "%s %s %s", sam_name, sam_rname, sam_cigar);
					int chr = ga.get_chromosome(string(sam_rname + 3));
					ga.set_chromosome_offset(chr, genome_len);
					genome_len += atoi(sam_cigar + 3);
				}
				continue;
			}
			else if (prev_name == "") {
				covx = new unsigned short[genome_len + 1];
				memset(covx, 0, sizeof(unsigned short) * (genome_len + 1));

				multiregion = new char[genome_len + 1];
				memset(multiregion, 0, (genome_len + 1));

				E("\t%5s %15s %15s %15s\n", "%%", "Partials", "Reads", "SAM lines");
			}
			sscanf(buffer, "%s %u %s %u %u %s %s %d %d %s", 
				sam_name, &sam_flag, sam_rname, &sam_pos, &sam_mapq, sam_cigar, sam_rnext, &sam_pnext, &sam_tlen, sam_read);
			sam_pos--;
			sam_pnext--;
		}

		if (prev_name != string(sam_name) || feof(fi)) {
			// iterate through all reads and obtain all pairs of partial transcripts
			set<read_pt_entry> result;
			
			#ifdef LOGIFY
			_sline1 = prev_name;
			_sline2 = prev_name + "/REV";
			#endif
		
			foreach (i, idx) 
				parse_read(i->first, i->second, result);
			#ifndef LOGIFY_ONLY
			// process ONLY if we have multi-mappings!
			if (result.size() > 1) { 
				if (read_id >= reads.size())
					reads.resize(read_id + 10000);
				foreach (i, result) {
					reads[read_id].entries.push_back(read::read_entry(
						PT(i->partial),
						make_pair(i->line1, i->line2),
						make_pair(i->start1, i->start2)
					));

					foreach (x, i->part1) 
				 		for (uint32_t _q=x->first;_q<= x->second;_q++)	multiregion[ga.get_chromosome_offset(i->chr1)+_q]=1;
				 	foreach (x, i->part2) 
				 		for (uint32_t _q=x->first;_q<= x->second;_q++)	multiregion[ga.get_chromosome_offset(i->chr2)+_q]=1;
				}
				read_id++;
			}
			// otherwise, just update the single-mapping partial counter
			else if (result.size() == 1) {
				//  ENSG00000135535.e6
				foreach (x, result.begin()->part1) 
				 	increase_coverage(result.begin()->chr1, x->first, x->second);
				foreach (x, result.begin()->part2) 
				 	increase_coverage(result.begin()->chr2, x->first, x->second);

				if (result.begin()->line1 != -1)
					single_maps.push_back(make_pair(result.begin()->line1, 
							result.begin()->partial->get_gene()));
				if (result.begin()->line2 != -1)
					single_maps.push_back(make_pair(result.begin()->line2, 
							result.begin()->partial->get_gene()));
			}
			// if there are no valid PTs, discard
			else {
				foreach (i, idx) { // crappy?! add ALL!
					if (i->second.line1 >= 0) crappy_reads.push_back(i->second.line1);
					if (i->second.line2 >= 0) crappy_reads.push_back(i->second.line2);
				}

				// coverage
				// only if not multi-crappy read
				if (idx.size() == 1) foreach (y, idx) {
					foreach (x, y->second.part1) 
						increase_coverage(y->first.chr1, x->first, x->second);
					foreach (x, y->second.part2) 
						increase_coverage(y->first.chr2, x->first, x->second);
				}
			}
			#endif

			prev_name = string(sam_name);
			idx.clear();
		}
		if (feof(fi)) break;
		// fix chromosome value
		if (string(sam_rnext) == "=")
			strcpy(sam_rnext, sam_rname);
		// get chromosome values
		int chr1 = ga.get_chromosome(sam_rname),
			chr2 = ga.get_chromosome(sam_rnext);
		// check chromosomes
		if (chr1 == -1 || chr2 == -1) { 
			crappy_reads.push_back(line);
			
			line++; continue; 
		}

		// is it first mate?
		if ((sam_flag & 0x8) || (sam_flag & 0x40)) {
			read_entry_key k(chr1, sam_pos,   
							 chr2, sam_pnext,    
							 abs(sam_tlen), sam_flag & 0x100); // not primary alignment
			multimap<read_entry_key, read_entry_value>::iterator i = idx.find(k);
			if (i != idx.end()) {
				i->second.line1 = line;
				parse_cigar(sam_pos, sam_cigar, sam_read, i->second.part1, i->second.indel1);
			}
			else {
				read_entry_value v;
				v.line1 = line;
				v.line2 = -1;
				parse_cigar(sam_pos, sam_cigar, sam_read, v.part1, v.indel1);
				idx.insert(make_pair(k, v));
			}
		}
		else {
			read_entry_key k(chr2, sam_pnext,
							 chr1, sam_pos,   
							 abs(sam_tlen), sam_flag & 0x100); 
			multimap<read_entry_key, read_entry_value>::iterator i = idx.find(k);
			if (i != idx.end()) {
				i->second.line2 = line;
				parse_cigar(sam_pos, sam_cigar, sam_read, i->second.part2, i->second.indel2);
			}
			else {
				read_entry_value v;
				v.line1 = -1;
				v.line2 = line;
				parse_cigar(sam_pos, sam_cigar, sam_read, v.part2, v.indel2);
				idx.insert(make_pair(k, v));
			}
		}

		read_length = max(read_length, (int)strlen(sam_read));
		line++;
		if (line % (1<<14) == 0) 
			E("\r\t%5.2lf %'15d %'15d %'15d", 
					100.0 * double(ftell(fi)) / f_size, PTs_set.size(), 
					reads.size(), line
			);
	}
	E("\n");

	sort(single_maps.begin(), single_maps.end());
	for(int i = 1;i < single_maps.size(); i++)
		assert(single_maps[i].first != single_maps[i-1].first);
	fclose(fi);

	delete[] buffer;
	delete[] sam_rname;
	delete[] sam_cigar;
	delete[] sam_name;
	delete[] sam_read;
	delete[] sam_rnext;
}
