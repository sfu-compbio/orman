/// 786

#ifndef PARTIAL_H__
#define PARTIAL_H__

#include "common.h"
#include "annotation.h"

using namespace std;

struct partial_transcript {
	genome_annotation::transcript *transcript;
	string signature;
	int    length;
	int 	 weight;

	partial_transcript (void) {}
	partial_transcript (genome_annotation::transcript *t, const string &s, int l, int w) : 
		transcript(t), signature(s), length(l), weight(w) {}

	bool operator== (const partial_transcript &s) const {
		return (transcript == s.transcript && signature == s.signature);
	}
	bool operator< (const partial_transcript &s) const {
		return (transcript < s.transcript || (transcript == s.transcript && signature < s.signature));
	}
};

struct read {
	struct read_key { // same line, it can map to multiple partial transcripts!
		const partial_transcript *partial;
		int line;

		read_key (void) {}
		read_key (int l, const partial_transcript *p) :
			partial(p), line(l) {}

		bool operator< (const read_key &r) const {
			return line < r.line || (line == r.line && partial < r.partial);
		}
	};
	struct read_entry {
		int chromosome;
		uint32_t position;
		const partial_transcript *partial;
		int partial_start;

		read_entry (void) {}
		read_entry (int c, uint32_t po, const partial_transcript *p, int s) : 
			chromosome(c),	position(po), partial(p), partial_start(s) {}
	};
	map<read_key, read_entry> entries; 
};

#endif // COMMON_H__

