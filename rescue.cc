/// 786

#include "common.h"
#include "interval.h"

using namespace std;

struct exon {
	uint32_t transcript_id,
				exon_id;
	uint32_t start, end;
	int 		chromosome;
	exon 		*prev, *next;

	exon() {}
	exon(uint32_t t, uint32_t e) :
		transcript_id(t), exon_id(e), next(0), prev(0) {}
	exon(uint32_t t, uint32_t e, uint32_t S, uint32_t E, int chr) :
		transcript_id(t), exon_id(e), start(S), end(E), next(0), chromosome(chr), prev(0) {}

	bool operator== (const exon &e) const {
		return (transcript_id == e.transcript_id && exon_id == e.exon_id);
	}
};
vector<exon> 			exons;
interval_tree<exon*> exon_tree_5p; // exon_tree_3p;
map<int, string>    	transcript_reverse_index;
map<string, int> 		chromosomes;

struct gene {
	uint32_t start, end;
	uint32_t unique;
	int chromosome;

	gene (void) : 
		start(-1), end(0), unique(0), chromosome(1) {}
	double expression (void) const {
		return double(unique) / (end - start + 1);
	}
};
map<string, gene> genes;
interval_tree<gene*> gene_tree;

/***********************************************************************/

void find_exons (const interval &iv, uint32_t flag, map<int, vector<exon*> > &result) {
	vector<exon*> olp;
	//! if (!(flag & 0x10))
	::exon_tree_5p.enumerate(iv, olp);
	//! else
	//!	::exon_tree_3p.enumerate(iv, olp);
	foreach (i, olp) 
		result[(*i)->transcript_id].push_back(*i);
}

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

void parse_gtf (const char *gtf_file) {
	FILE * fi = fopen(gtf_file, "r");

	char 		buffer[MAX_BUFFER];
	char 		gtf_seqname[MAX_BUFFER],
		  		gtf_source[MAX_BUFFER],
		  		gtf_feature[MAX_BUFFER];
	uint32_t gtf_start, gtf_end;
	char 		gtf_score, gtf_strand, gtf_frame;

	map<string, int> transcripts;
	map<string, int>::iterator ti;
	int transcripts_count = 0;
	map<string, string> attributes;

	while (fgets(buffer, MAX_BUFFER, fi)) {
		int of; sscanf(buffer, "%s %s %s %u %u %c %c %c%n",
				gtf_seqname, gtf_source, gtf_feature,
				&gtf_start, &gtf_end, &gtf_score, &gtf_strand, &gtf_frame,
				&of);

		if (!strcmp(gtf_feature, "exon")) {
			auto it = chromosomes.find(gtf_seqname);
			int chromosome;
			if (it != chromosomes.end()) 
				chromosome = it->second;
			else 
				chromosomes.insert(make_pair(gtf_seqname, chromosome = chromosomes.size()));

			gene &g = genes[attributes["gene_id"]];
			g.start = min(g.start, gtf_start);
			g.end   = max(g.end,   gtf_end);
			g.chromosome = chromosome;

			parse_gtf_attributes(buffer + of, attributes);
			string transcript_uid = attributes["gene_id"] + "+" + attributes["transcript_id"];

			int transcript_id;
			if ((ti = transcripts.find(transcript_uid)) != transcripts.end())
				transcript_id = ti->second;
			else {
				transcript_id = transcripts_count;
				transcripts[transcript_uid] = transcripts_count;
				transcript_reverse_index[transcripts_count++] = transcript_uid;
			}
			int exon_id = atoi(attributes["exon_number"].c_str());

			/*E("Found exon\n\t%s\n\t%u->%u, G %s, T %s, T# %d, E %s->%d   %c\n",
				buffer,
				gtf_start, gtf_end, attributes["gene_id"].c_str(),
				attributes["transcript_id"].c_str(), transcript_id,
				attributes["exon_number"].c_str(), exon_id, gtf_strand);
			W();*/

			//! if (gtf_strand == '+')
			exons.push_back(exon(transcript_id, exon_id, gtf_start, gtf_end, chromosome));
		}
	}
	for (int i = 0; i < exons.size(); i++) {
		if (i) {
			exons[i - 1].next = &exons[i];
			exons[i].prev = &exons[i - 1];
		}
		exon_tree_5p.insert(interval(exons[i].start, exons[i].end), &exons[i]);
	}
	foreach (g, genes)
		gene_tree.insert(interval(g->second.start, g->second.end), & g->second);

	fclose(fi);
}

gene * getgene(int i, int chr) {
	vector<gene*> gx;
	gene_tree.enumerate(i, gx);
	
	foreach(g, gx) {
		if ((*g)->chromosome == chr)
			return (*g);
	}
	return 0;
//	if(gx.size()>1)E("UUUUPS %d MULTIGENE MWHAHAHAHAH!\n", i);
//	return gx[0];
}

bool parse_cigar (uint32_t start_pos, const char *cigar, vector<interval> &result) {
	int num = 0;
	while (*cigar) {
		if (isdigit(*cigar))
			num = 10 * num + (*cigar - '0');
		else if (*cigar == 'M') {
			result.push_back(interval(start_pos, start_pos + num));
			start_pos += num;
			if (num < 5)
				return 0;
			num = 0;
		}
		else if (*cigar == 'N' || *cigar == 'D') {
			start_pos += num;
			num = 0;
		}
		else if (*cigar == 'I' || *cigar == 'S') {
		//	start_pos -= num;
			num = 0;
		}
		// H, O, =, X
		else {
			E("CIGAR feature %c not implemented!", *cigar);
			return 0;
		}
		cigar++;
	}
	return 1;
}

/*******************************************************************************/

struct semitranscript {
	int		transcript_id;
	string	signature;
	int      length;
	int 		weight;

	semitranscript () {}
	semitranscript (int t, const string &s) : 
		transcript_id(t), signature(s), length(0), weight(1) {}

	bool operator< (const semitranscript &s) const {
		return (transcript_id < s.transcript_id || (transcript_id == s.transcript_id && signature < s.signature));
	}
};
map<semitranscript, pair<int, int> > sets;
int sets_count = 0;

string numtostr (int n) {
	string s = "";
	while (n) {
		s = char(n % 10 + '0') + s;
		n /= 10;
	}
	return s;
}

int get_semitranscript (const char *read, int transcript_id, const vector<exon*> &exons, const vector<interval> &parts, int &st_pos) {
	
	semitranscript st(transcript_id, "");
	int pi = 0, pex = -1;
	bool prev_intron = false;

	st_pos = -1;
	bool first = true;
	if (parts[0].first >= exons[0]->start)
		st_pos = parts[0].first - exons[0]->start;
	else 
		st_pos = -1;

	st.length = 0; 
	bool is_intronic = 0;
	foreach (ep, exons) {
		exon *e = *ep;
		if (!prev_intron && parts[pi].first < e->start && !(pi == 0 && e->start - parts[pi].first < 5))  {
			is_intronic = 1;
		}
		if (pex == e->exon_id) {
			return -1;
		}
		pex = e->exon_id;
		while (pi < parts.size() && parts[pi].second <= e->end)
			pi++;
		if (prev_intron = (pi < parts.size () && parts[pi].first < e->end && !(pi == parts.size() - 1 && parts[pi].second - e->end < 5))) {
			is_intronic = 1;
		}
	}
	if (is_intronic)  // forget them for now
		return -1;

	st.weight = 1;
	exon *e = exons[0];
	int dir = (exons.size() > 1 && exons[0]->exon_id > exons[1]->exon_id);
	int skips = 0;
	for (int i = 1; i < exons.size(); i++) {
		e = dir ? e->prev : e->next;
		if (e != exons[i]) {
		//	E("(%d) ",e->exon_id);
			skips++;
			while (e != exons[i]) {
				st.weight += e->end - e->start + 1;
				e = dir ? e->prev : e->next;
			}
		}
		// E("%d\n", exons[i]->exon_id);
	}
	if (skips >= 2) {
		E("Bigskip! %d %s\n", st.transcript_id, st.signature.c_str());
	}

/*	
	typeof(sets.begin()) it = sets.find(st);
	if (it != sets.end()) {
		it->second.second++;
		return it->second.first;
	}
	else {
		sets[st] = make_pair(sets_count++, 1);
		return sets_count - 1;
	}*/
	return 10;
}

/****************************************************************************/

void parse_sam (const char *sam_file) {
	FILE *fi = fopen(sam_file, "r");

	char 		buffer[MAX_BUFFER];
	char 		sam_name[MAX_BUFFER],
		  		sam_rname[MAX_BUFFER],
		  		sam_cigar[MAX_BUFFER],
				sam_read[MAX_BUFFER],
				sam_rnext[MAX_BUFFER];
	uint32_t sam_flag, sam_pos;
	int32_t  sam_pnext, sam_tlen;
	uint8_t  sam_mapq;

	// read single map
	E("Reading single mapped reads ...\n");
	while (fgets(buffer, MAX_BUFFER, fi)) {
		if (buffer[0] == '@') continue;

		char *cx;
		if ((cx = strstr(buffer, "NH:i:1")) != 0 && !isdigit(cx[6])) {
			sscanf(buffer, "%s %u %s %u %u %s %s %d %d %s",
				sam_name, &sam_flag, sam_rname,
				&sam_pos, &sam_mapq, sam_cigar, 
				sam_rnext, &sam_pnext, &sam_tlen, sam_read);

			auto chi = chromosomes.find(sam_rname);
			if (chi == chromosomes.end()) {
				E("Chromosome %s not found in GTF!\n", sam_rname);
				continue;
			}

			gene *g = getgene(sam_pos, chi->second);
			if (g) g->unique++;
			continue;
		}
	}

	E("Reading multimap reads...\n");
	rewind(fi);
	map<string, pair<int, int> > read_genes;
	while (fgets(buffer, MAX_BUFFER, fi)) {
		if (buffer[0] == '@') continue;

		char *cx;
		if ((cx = strstr(buffer, "NH:i:1")) != 0 && !isdigit(cx[6])) 
			continue;
		
		sscanf(buffer, "%s %u %s %u %u %s %s %d %d %s",
			sam_name, &sam_flag, sam_rname,
			&sam_pos, &sam_mapq, sam_cigar, 
			sam_rnext, &sam_pnext, &sam_tlen, sam_read);

		int l = strlen(sam_name);
		sam_name[l++] = '/';
		sam_name[l++] = ((sam_flag & 0x40) ? '1' : '2');
		sam_name[l] = 0;

		// single?!
		auto chi = chromosomes.find(sam_rname);
		if (chi == chromosomes.end()) {
			E("Chromosome %s not found in GTF!\n", sam_rname);
			continue;
		}
		vector<interval> parts;
		if (!parse_cigar(sam_pos, sam_cigar, parts)) 
			continue;
		map<int, vector<exon*> > transcripts;
		foreach (i, parts)  
			find_exons(*i, sam_flag, transcripts);	
		
		char ok = 0;
		foreach (t, transcripts) 
			if (t->second.size() && t->second[0]->chromosome == chi->second) {
				int st_pos;
				int s_id = get_semitranscript(sam_name, t->first, t->second, parts, st_pos);
				ok |= (s_id != -1);
			}	
		if (!ok) continue;

		/************************/

		gene *g = getgene(sam_pos, chi->second);
		if (!g) continue;
		auto git = read_genes.find(sam_name);
		if (git == read_genes.end() ||
				getgene(sam_pos, chi->second)->expression() > getgene(git->second.first, git->second.second)->expression())
			read_genes[sam_name] = make_pair(sam_pos, chi->second);
	}

	rewind(fi);

	int single = 0,
		 destroyed = 0, 
		 discarded  = 0,
		 total = 0;

	while (fgets(buffer, MAX_BUFFER, fi)) {
		if (buffer[0] == '@')  {
			fputs(buffer, stdout);
			continue;
		}
		total++;
		char *cx;
		if ((cx = strstr(buffer, "NH:i:1")) != 0 && !isdigit(cx[6])) { 
			fputs(buffer, stdout);
			single++;
			continue;
		}

		sscanf(buffer, "%s %u %s %u %u %s %s %d %d %s",
				sam_name, &sam_flag, sam_rname,
				&sam_pos, &sam_mapq, sam_cigar, 
				sam_rnext, &sam_pnext, &sam_tlen, sam_read);
		int l = strlen(sam_name);
		sam_name[l++] = '/';
		sam_name[l++] = ((sam_flag & 0x40) ? '1' : '2');
		sam_name[l] = 0;

		auto itx = read_genes.find(sam_name);
		if (itx == read_genes.end()) {
			destroyed++;
			continue;
		}
		auto &x = itx->second;
		if (x.second > -1 && sam_pos == x.first) {
			fputs(buffer, stdout);
			x.second = -1;
		}
		else discarded++;
	}

	E("Read %d single %d destroyed %d discarded %d\n", total, single, destroyed, discarded);

	fclose(fi);
}

int main (int argc, char **argv) {
	setlocale(LC_ALL, "");
	E("Behold! THE HEPEK is starting! Usage ! [gtf] [sam]\n");

	char buffer[MAX_BUFFER];
	zaman_last();
	E("Parsing GTF file %s ... ", realpath(argv[1], buffer));
	parse_gtf(argv[1]);
	E("done in %d seconds!\n", zaman_last());
	E("Parsing SAM file %s ... ", realpath(argv[2], buffer));
	parse_sam(argv[2]);
	E("done in %d seconds!\n", zaman_last());

	return 0;
}

