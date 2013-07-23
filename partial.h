/// 786

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

	PT_single (void) {}
	PT_single (genome_annotation::transcript *t, const string &s, int l, int w) : 
		transcript(t), signature(s), length(l), weight(w) {}

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

struct PT {
	PTsp first;
	PTsp second;

	PT() : first(0), second(0) {}
	PT(PTsp a, PTsp b): first(a), second(b) {}
	bool operator== (const PT &x) const { return first==x.first && second==x.second; }
	bool operator< (const PT &x) const { return first<x.first || (first==x.first && second<x.second); }
};
//ypedef pair<PTsp, PTsp> 	PT;  
int get_single_coverage(const PT &p);
int get_single_coverage(const PT &p, int k);

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
//	int solution;

//	read () : solution(0) {}
};

#endif // COMMON_H__

