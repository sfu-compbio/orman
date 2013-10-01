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

#include "annotation.h"
using namespace std;

inline const char *parse_empty (const char *c) {
	while (*c && isspace(*c)) c++;
	return c;
}

inline const char *parse_field (const char *c, char *d) {
	while (*c && !isspace(*c)) *d++ = *c++;
	if (*(d - 1) == ';' && *(d - 2) == '"') d = d - 2; // TODO: DANGEROUS!
	*d = 0;
	return c;
}

void parse_gtf_attributes (const char *attributes, map<string, string> &result) {
	char name[MAX_BUFFER], value[MAX_BUFFER];
	while (attributes) {
		attributes = parse_empty(attributes);
		attributes = parse_field(attributes, name);
		if (*name == 0) break;
		attributes = parse_empty(attributes);
		attributes = parse_field(attributes, value);
		result[string(name)] = string(value + (value[0] == '"'));
	}
}

template<typename T>
T &get_or_insert(map<string, T> &m, const string &k) {
	T &t = m[k];
	return t;
}

void genome_annotation::parse_gtf (const char *gtf_file) {
	FILE * fi = fopen(gtf_file, "r");

	char 		buffer[MAX_BUFFER];
	char 		gtf_seqname[MAX_BUFFER],
		  		gtf_source[MAX_BUFFER],
		  		gtf_feature[MAX_BUFFER];
	uint32_t 	gtf_start, gtf_end;
	char 		gtf_score, gtf_strand, gtf_frame;

	
	map<string, string> attributes;
	while (fgets(buffer, MAX_BUFFER, fi)) {
		int of; sscanf(buffer, "%s %s %s %u %u %c %c %c%n",
				gtf_seqname, gtf_source, gtf_feature,
				&gtf_start, &gtf_end, &gtf_score, &gtf_strand, &gtf_frame,
				&of);
		if (!strcmp(gtf_feature, "exon")) {
			map<string, int>::iterator it = chromosomes.find(gtf_seqname);
			int chromosome;
			if (it != chromosomes.end()) 
				chromosome = it->second;
			else { 
				chromosome = chromosomes.size();
				chromosomes.insert(make_pair(gtf_seqname, chromosome));
			}

			parse_gtf_attributes(buffer + of, attributes);

			gene &g = get_or_insert(genes, attributes["gene_id"]);
			transcript &t = get_or_insert(g.transcripts, attributes["gene_id"]);

			if (attributes.find("partial_ex") == attributes.end())
				throw string("GTF file does not contain partial_ex field. In order to use this GTF file, you need to process it with the provided ormanGTF script.");

			// Make it 0-based
			gtf_start--;
			gtf_end--; // Interval right is INCLUSIVE here!
			
			t.exons.push_back(exon(gtf_start, gtf_end, attributes["partial_ex"]));
			g.chromosome = chromosome;
			g.name = attributes["gene_id"];
		}
	}
	transcript_count = 0;
	foreach (gi, genes) 
		foreach (ti, gi->second.transcripts) {
			transcript &t = ti->second;

			t.gene = &gi->second;
			t.name = ti->first;
			t.id   = transcript_count++;
			t.length_ = 0;

			sort(t.exons.begin(), t.exons.end());
			for (int i = 1; i < t.exons.size(); i++) 
				if (t.exons[i].sid == t.exons[i - 1].sid) {
					t.exons.erase(t.exons.begin() + i);
					i--;
				}
			for (int i = 0; i < t.exons.size(); i++) {
				t.exons[i].id = i;
				t.exons[i].transcript = &t;
				t.length_ += t.exons[i].end-t.exons[i].start+1;

				exon_tree.insert(interval(t.exons[i].start, t.exons[i].end), &t.exons[i]);
			}
		}
	chromosomes.insert(make_pair("*", chromosomes.size()));
	fclose(fi);
}

void genome_annotation::get_exons (const interval &i, vector<exon*> &exons) {
	exon_tree.enumerate(i, exons);
}

genome_annotation::gene *genome_annotation::find_gene (uint32_t pos, int chr) {
	vector<exon*> ex;
	exon_tree.enumerate(pos, ex);
	
	foreach (e, ex) {
		if ((*e)->transcript->gene->chromosome == chr) 
			return (*e)->transcript->gene;
	}

	ex.clear();
	exon_tree.enumerate(pos - 4, ex);
	exon_tree.enumerate(pos + 4, ex);
	foreach (e, ex) {
		if ((*e)->transcript->gene->chromosome == chr)
			return (*e)->transcript->gene;
	}

	return 0;
}

