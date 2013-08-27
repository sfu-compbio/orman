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

/*struct PT {
	PTsp first;
	PTsp second;

	PT() : first(0), second(0) {}
	PT(PTsp a, PTsp b): first(a), second(b) {}
	bool operator== (const PT &x) const { return first==x.first && second==x.second; }
	bool operator< (const PT &x) const { return first<x.first || (first==x.first && second<x.second); }
};*/
//ypedef pair<PTsp, PTsp> 	PT;  
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
//	int solution;

//	read () : solution(0) {}
};

#endif // COMMON_H__

