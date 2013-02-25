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

/*struct partial_transcript {
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
};*/

typedef pair<const partial_transcript_single*, const partial_transcript_single*> partial_transcript;  
struct read {
	struct read_entry {
		partial_transcript partial;
		pair<uint32_t, uint32_t> position, partial_start;
		pair<int, int> line;
		
		read_entry (void) {}
		read_entry (const partial_transcript &p, const pair<int, int> &l, const pair<uint32_t, uint32_t> &po, const pair<uint32_t, uint32_t> &s) : 
			partial(p), line(l), position(po), partial_start(s) {}
	};
	vector<read_entry> entries;
};

#endif // COMMON_H__

