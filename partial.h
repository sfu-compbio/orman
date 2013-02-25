/// 786

#ifndef PARTIAL_H__
#define PARTIAL_H__

#include "common.h"
#include "annotation.h"

using namespace std;

struct partial_transcript_single {
	genome_annotation::transcript *transcript;
	string signature;
	int    length;
	int 	 weight;

	partial_transcript_single (void) {}
	partial_transcript_single (genome_annotation::transcript *t, const string &s, int l, int w) : 
		transcript(t), signature(s), length(l), weight(w) {}

	bool operator== (const partial_transcript_single &s) const {
		return (transcript == s.transcript && signature == s.signature);
	}
	bool operator< (const partial_transcript_single &s) const {
		return (transcript < s.transcript || (transcript == s.transcript && signature < s.signature));
	}
};

struct partial_transcript {
	const partial_transcript_single *p1, *p2;

	partial_transcript(const partial_transcript_single *x1, const partial_transcript_single *x2) :
		p1(x1), p2(x2) {}

	bool operator< (const partial_transcript &s) const {
		return (!p2 || p1<s.p1 || (p1==s.p1&&p2<s.p2));
	}

	int length () const {
		return p1->length;
	}
	int weight () const {
		if (!p2)
			return 1000000 + p1->weight;
		else if (p1->transcript == p2->transcript)
			return p1->weight + p2->weight;
		else 
			return p1->weight + p2->weight + 100000;
	}
};

struct read {
	struct read_key { // same line, it can map to multiple partial transcripts!
		const partial_transcript *partial;
		int line, line2;

		read_key (void) {}
		read_key (int l, int l2, const partial_transcript *p) :
			partial(p), line(l), line2(l2) {}

		bool operator< (const read_key &r) const {
			return line < r.line || (line == r.line && partial < r.partial);
		}
	};
	struct read_entry {
		uint32_t position;
		uint32_t partial_start;

		read_entry (void) {}
		read_entry (uint32_t po, uint32_t s) : 
			position(po), partial_start(s) {}
	};
	map<read_key, pair<read_entry, read_entry> > entries; 
};

#endif // COMMON_H__

