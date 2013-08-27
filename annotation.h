/// 786

#ifndef ANNOTATION_H__
#define ANNOTATION_H__

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
			//return exons[exons.size() - 1].end - exons[0].start + 1;
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
		auto it = chromosomes.find(s);
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


