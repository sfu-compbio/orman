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
			start(s), end(e), sid(i) {}
		bool operator< (const exon &e) const {
			return start < e.start;
		}
	};
	struct transcript {
		int id;
		string name;
		vector<exon> exons;
		struct gene *gene;

		transcript (void) {}

		uint32_t length (void) const {
			return exons[exons.size() - 1].end - exons[0].start + 1;
		}

		int chromosome (void) const {
			return gene->chromosome;
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
	map<string, int> chromosomes;

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


