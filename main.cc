/// 786

#include "common.h"
#include "annotation.h"
#include "interval.h"
#include "rescue.h"
#include "partial.h"
#include "orman.h"
#include <tuple>
#include <cstdarg>
using namespace std;

#ifdef LOGIFY_ONLY
#define LOGIFY
#endif

const char *DEBUG_LOG_FILE 	= "run/orman.dbg";
const char *LOG_FILE = "run/orman.log";
char junk_file_path[MAX_BUFFER];
bool retain_junk = false;
int read_length = 0;
const char *gene_sam_flag = "YG";
const char *partial_sam_flag = "YP";

FILE *flog;
#define LOG(c,...) fprintf(flog,c,##__VA_ARGS__)

const char *S (const char* f, ...) {
	static char bf[MAX_BUFFER];
	va_list args;
	va_start(args, f);
	vsprintf(bf, f, args);
	va_end(args);
	return (bf);
}

typedef genome_annotation::transcript transcript;
typedef genome_annotation::exon exon;

struct indel {
	enum indel_type { INSERT, DELETE };

	int first, second;
	indel_type type;
	string insert;

	indel (void) {}
	indel (int f, int s, indel_type t) :
		first(f), second(s), type(t) {}
};

/*******************************************************************************/

set<PTs> PTs_set;
vector<struct read> reads;
genome_annotation ga;

vector< vector<short> > covx;

int get_single_coverage(const PT &p, int k) {
	return covx[p.first->transcript->chromosome()]
			   [p.first->transcript->position(k)];
}

void increase_coverage(int chr, int p1, int p2) {
	assert(chr < covx.size());
	assert(p1 <= p2);
	assert(p2 < covx[chr].size());
	for (int i = p1; i < p2; i++)
		covx[chr][i]++;	
}

/*******************************************************************************/

bool parse_cigar (uint32_t start_pos, const char *cigar, const char *read, vector<interval> &result, vector<indel> &indels) {
	int num = 0;

	if (strcmp(cigar, "*") == 0)
		return 0;
	
	while (*cigar) {
		if (isdigit(*cigar))
			num = 10 * num + (*cigar - '0');
		else {
			if (*cigar == 'M') {
				result.push_back(interval(start_pos, start_pos + num));
				start_pos += num;
				read += num;
			}
			else if (*cigar == 'D') { // delete FROM REFERENCE GENOME
				indels.push_back(indel(start_pos, start_pos + num, indel::indel_type::DELETE));
			//	start_pos += num;
			}
			else if (*cigar == 'N') // TODO check is the same as DEL
				start_pos += num;
			else if (*cigar == 'I') { // insert IN READ 
				indel x(start_pos, start_pos + num, indel::indel_type::INSERT);
				x.insert = string(read, num);
				indels.push_back(x);
				read += num;
			}
			else if (*cigar == 'S') {
				read += num;
			}
			else { // H O = X
				E("CIGAR feature %c not implemented!", *cigar);
				return 0;
			}
			num = 0;
		}
		cigar++;
	}
	return 1;
}

/*******************************************************************************/

#ifdef LOGIFY
	static string get_partial_transcript_error;
#endif

PTsp get_partial_transcript (transcript *t, const vector<exon*> &exons, const vector<interval> &parts, const vector<indel> &indels, int &starting_position) {
	// disregard junk
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

	// indels, exons are sorted; detect them
	string signature = "";
	int length = 0;
	int weight = 1;

	int indel_idx = 0;
	for (int i = 0; i < exons.size(); i++) {
		exon *e = exons[i];
		// penalize missing exons by their length
		if (i) for (int j = exons[i - 1]->id + 1; j < exons[i]->id; j++)
			weight += 100 * (t->exons[j].end - t->exons[j].start + 1);
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
			if (indels[indel_idx].type == indel::indel_type::INSERT) { // reference insert
				signature += "I." + indels[indel_idx].insert;
				length += indels[indel_idx].second - indels[indel_idx].first + 1;
			}
			else { // reference delete
				signature += "D";
				length += indels[indel_idx].second - indels[indel_idx].first + 1;
			}
			signature += "," + numtostr(indels[indel_idx].first - e->start) // start
						  + "," + numtostr(indels[indel_idx].second - indels[indel_idx].first) // length
						  + "]";
			// penalize indels
			weight += 10000;
			indel_idx++;
		}
	}

	PTs pt(t, signature, length, weight);
	auto pti = PTs_set.insert(pt);
	return &(* pti.first);
}

/****************************************/

void trimx (vector<interval> &parts) {
	if (parts.size() == 0)
		return;
	if (parts[0].second - parts[0].first + 1 < 4)
		parts = vector<interval>(parts.begin() + 1, parts.end());
	if (parts[parts.size() - 1].second - parts[parts.size() - 1].first + 1 < 4)
		parts.pop_back();
}

void makecand (vector<interval> parts, int chromosome, map<transcript*, vector<exon*> > &candidates) {
	if (chromosome == -1) return;
	trimx(parts);

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

#ifdef LOGIFY
string _sline2, _sline1;
#endif

/*******************************************************************************/

vector<pair<int64_t, genome_annotation::gene*>> single_maps;
vector<int64_t> crappy_reads;

struct read_entry_key {
	int32_t  chr1, chr2;
	uint32_t pos1, pos2;

	int32_t  tlen;
	char     is_secondary;

	read_entry_key() {}
	read_entry_key(int32_t c1, uint32_t p1, int32_t c2, uint32_t p2, int32_t t, char i):
		chr1(c1), pos1(p1), chr2(c2), pos2(p2), tlen(t), is_secondary(i) {}
	bool operator< (const read_entry_key& x) const {
		if (chr1 != x.chr1) return chr1 < x.chr1;
		if (chr2 != x.chr2) return chr2 < x.chr2;
		if (pos1 != x.pos1) return pos1 < x.pos1;
		if (pos2 != x.pos2) return pos2 < x.pos2;
		if (tlen != x.tlen) return tlen < x.tlen;
		return is_secondary < x.is_secondary;
	}
};

struct read_entry_value {
	int64_t line1, line2;
	vector<interval> part1, part2;
	vector<indel> indel1, indel2;
};

struct read_pt_entry {
	PTsp partial1, partial2;
	int chr1, chr2;
	int64_t line1, line2;
	int32_t start1, start2;
	vector<interval> part1, part2; 

	read_pt_entry() {}
	read_pt_entry(PTsp p1, int c1, int64_t l1, int32_t s1, const vector<interval> &i1,
				  PTsp p2, int c2, int64_t l2, int32_t s2, const vector<interval> &i2):
		partial1(p1), chr1(c1), line1(l1), start1(s1), part1(i1),
		partial2(p2), chr2(c2), line2(l2), start2(s2), part2(i2) {}
	bool operator< (const read_pt_entry& x) const {
		return (line1 < x.line1) || (line1 == x.line1 && line2 < x.line2);
	}
};

/*
 * get all pairs of partial transcripts for each paired read
 */
void parse_read (const read_entry_key &rk, const read_entry_value &rv, set<read_pt_entry> &result) {
	map<transcript*, vector<exon*> > candidates1, candidates2;
	makecand(rv.part1, rk.chr1, candidates1);
	makecand(rv.part2, rk.chr2, candidates2);

	//// crappy coverage
	// #ifdef LOGIFY
	// LOG("%s ", _sline1.c_str());
	// foreach (ci1, candidates1) {
	// 	int starting_position1 = 0; 
	// 	PTsp pt1 = get_partial_transcript(ci1->first, ci1->second, parts1, indels1, starting_position1);
	// 	if (pt1) LOG("%s ", string(pt1->transcript->gene->name /*+ "." + pt1->transcript->name */+ ":" + pt1->signature).c_str());
	// 	else LOG("%s ", get_partial_transcript_error.c_str());
	// }
	// LOG("\n");
	// LOG("%s ", _sline2.c_str());
	// foreach (ci2, candidates2) {
	// 	int starting_position2 = 0; 
	// 	PTsp pt2 = get_partial_transcript(ci2->first, ci2->second, parts2, indels2, starting_position2);
	// 	if (pt2) LOG("%s ", string(pt2->transcript->gene->name + /*"." + pt2->transcript->name + */":" + pt2->signature).c_str());
	// 	else LOG("%s ", get_partial_transcript_error.c_str());
	// }
	// LOG("\n");
	// #endif

	bool crappy = true;
	foreach (ci1, candidates1) {
		int start1 = 0; 
		PTsp pt1 = get_partial_transcript(ci1->first, ci1->second, rv.part1, rv.indel1, start1);
		// first mate HAS TO BE VALID!
		if (pt1) {
			crappy = false;
			// special case if there is no valid pair
			bool no_null = false;
			if (candidates2.size()) foreach (ci2, candidates2) {
				int start2 = 0;
				PTsp pt2 = get_partial_transcript(ci2->first, ci2->second, rv.part2, rv.indel2, start2);
				if (pt2) result.insert(read_pt_entry(pt1, rk.chr1, rv.line1, start1, rv.part1,
													 pt2, rk.chr2, rv.line2, start2, rv.part2));
				no_null |= (pt2 != 0);
			}
			if (!no_null) result.insert(read_pt_entry(pt1, rk.chr1, rv.line1, start1, rv.part1,
													    0, rk.chr2, rv.line2,      0, rv.part2));
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

	char buffer[MAX_BUFFER];
	char sam_name[MAX_BUFFER],
		 sam_rname[MAX_BUFFER],
		 sam_cigar[MAX_BUFFER],
		 sam_read[MAX_BUFFER],
		 sam_rnext[MAX_BUFFER];
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
					int len = atoi(sam_cigar + 3);
					if (chr > covx.size()) covx.resize(chr + 10);
					covx[chr].resize(len, 0);
					E("\tChromosome %2s [%02d] of size %'12d\n", sam_rname+3, chr, len);
				}
				continue;
			}
			else if (prev_name == "") {
				E("\t%5s %15s %15s %15s\n", "%%", "Partials", "Reads", "SAM lines");
			}
			sscanf(buffer, "%s %u %s %u %u %s %s %d %d %s", 
				sam_name, &sam_flag, sam_rname, &sam_pos, &sam_mapq, sam_cigar, sam_rnext, &sam_pnext, &sam_tlen, sam_read);
		}

		// E(">%d\n",line);
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
				foreach (i, result) 
					reads[read_id].entries.push_back(read::read_entry(
						PT(i->partial1, i->partial2),
						make_pair(i->line1, i->line2),
						make_pair(i->start1, i->start2)
					));
				read_id++;
			}
			// otherwise, just update the single-mapping partial counter
			else if (result.size() == 1) {
				// coverage!
				foreach (x, result.begin()->part1) 
				 	increase_coverage(result.begin()->chr1, x->first, x->second);
				foreach (x, result.begin()->part2) 
				 	increase_coverage(result.begin()->chr2, x->first, x->second);
			
				if (result.begin()->line1 != -1)
					single_maps.push_back(make_pair(result.begin()->line1, 
							result.begin()->partial1->get_gene()));
				if (result.begin()->line2 != -1)
					single_maps.push_back(make_pair(result.begin()->line2, 
							result.begin()->partial2 ? result.begin()->partial2->get_gene() : 0));
			}
			// if there are no valid PTs, discard
			else {
				if (retain_junk) foreach (i, idx) { // crappy?! add ALL!
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
			auto i = idx.find(k);
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
			auto i = idx.find(k);
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
	for(int i=1;i<single_maps.size();i++) {
		if(single_maps[i].first==single_maps[i-1].first) {
			E("OOOOOOOOOOOOOps! %d %s %d %s\n", 
				single_maps[i].first,single_maps[i].second?single_maps[i].second->name.c_str():"_",
				single_maps[i-1].first,single_maps[i-1].second?single_maps[i-1].second->name.c_str():"_"
				);
			abort();
		}
	}

	fclose(fi);
}

void resolve () {
	foreach (r, reads) {
		if (r->entries.size() == 1) {
			assert(r->entries.begin()->line.first != -1);
			single_maps.push_back(make_pair(
						r->entries.begin()->line.first, 
						r->entries.begin()->partial.first->get_gene()));
			if (r->entries.begin()->line.second != -1)
				single_maps.push_back(make_pair(
						r->entries.begin()->line.second, 
						r->entries.begin()->partial.second ? r->entries.begin()->partial.second->get_gene() : 0));
		}
	}
	reads.clear();
	sort(single_maps.begin(), single_maps.end());
	if (retain_junk)
		sort(crappy_reads.begin(), crappy_reads.end());
	
	for(int i=1;i<single_maps.size();i++) {
		if(single_maps[i].first==single_maps[i-1].first)
		{
			E("OOOOOOOOOOOOOps!");
			exit(1);
		}
	}
}

void write_sam (const char *old_sam, const char *new_sam) {
	FILE *fi = fopen(old_sam, "r"),
		 *fo = fopen(new_sam, "w");

	FILE *fj;
	if (!retain_junk || !strcmp(junk_file_path, "="))
		fj = fo;
	else
		fj = fopen(junk_file_path, "w");

	// obtain the file size
	fseek(fi, 0, SEEK_END);
	int64_t f_size = ftell(fi);
	fseek(fi, 0, SEEK_SET);


	char buffer[MAX_BUFFER];
	int wr = 0,
		 di = 0,
		 i = 0,
		 line = 0;

	int j = 0;
	int ore = 0;

	char nam[MAX_BUFFER];
	string rev_name = "";

	E("\t%5s %15s\n", "%%", "SAM lines");
	while (fgets(buffer, MAX_BUFFER, fi)) {
		if (buffer[0] == '@') {	
			fputs(buffer, fo);
			continue;
		}
		sscanf(buffer, "%s", nam);
		if (string(nam) != rev_name) { ore++; rev_name = string(nam); }
		// TODO fix NH:i:...
		if (i < single_maps.size() && line == single_maps[i].first) {
			int l = strlen(buffer) - 1;
			sprintf(buffer + l, "\t%s:Z:%s\n", gene_sam_flag, 
					single_maps[i].second ? single_maps[i].second->name.c_str() : "_"
					);
			fputs(buffer, fo);
			i++;
			wr++;
		}
		else if (retain_junk && j < crappy_reads.size() && line == crappy_reads[j]) {
			int l = strlen(buffer) - 1;
			sprintf(buffer + l, "\t%s:Z:_\n", gene_sam_flag);
			fputs(buffer, fj);
			j++;
		}
		else {
			char name[MAX_BUFFER]; int sam_flag;
			sscanf(buffer, "%s %d", name, &sam_flag);
//			L("%s/%d at %d discarded\n", name, ((sam_flag & 0x8) || (sam_flag & 0x40)) ? "1" : "2", line);
			di++;
		}
		line++;
		if (line % (1<<14) == 0) 
			E("\r\t%5.2lf %'15d", 100.0 * double(ftell(fi)) / f_size, line);
	}
	E("\nOK!\twritten %'d, discarded %'d, total %'d\n", wr, di, line);
	E("\tcrappy  %'d, total %'d\n", j, wr + j);
	E("\tread count %'d\n", ore);

	fclose(fo);
	fclose(fi);
	if (fj != fo) fclose(fj);
}

/*******************************************************************************/

void parse_opt (int argc, char **argv, char *gtf, char *sam, char *newsam, char *mode) {
	int opt; 
	struct option long_opt[] = {
		{ "help",    0, NULL, 'h' },
		{ "gtf",  	 1, NULL, 'g' },
		{ "junk",  	 1, NULL, 'j' },
		{ "sam",	 1, NULL, 's' },
	//	{ "mode",	 1, NULL, 'm' },
		{ NULL,     0, NULL,  0  }
	};
	const char *short_opt = "hj:g:s:";
	do {
		opt = getopt_long (argc, argv, short_opt, long_opt, NULL);
		switch (opt) {
			case 'h':
				exit(0);
			case 'g':
				strncpy(gtf, optarg, MAX_BUFFER);
				break;
			case 'j':
				strncpy(junk_file_path, optarg, MAX_BUFFER);
				retain_junk = true;
				break;
			case 's':
				strncpy(sam, optarg, MAX_BUFFER);
				break;
		//	case 'm':
		//		strncpy(mode, optarg, MAX_BUFFER);
		//		break;
			case -1:
				break;
			default: {
				exit(1);
			}
		}
	} while (opt != -1);
	strncpy(newsam, argv[optind], MAX_BUFFER);
}

int main (int argc, char **argv) {
	setlocale(LC_ALL, "");
	char buffer[MAX_BUFFER];

	E("Behold! Former Uniqorn Republic of Orman is starting!\n\tUsage: orman -g[gtf] -s[sam] -m[mode:(single|rescue|orman)] [new_sam]\n");
	E("\tCompile time: %s\n", COMPILE_TIME);
	#ifdef LOGIFY
		E("\tLog status: enabled\n");
		flog = fopen(LOG_FILE, "w");
		E("\tLog file: %s\n", realpath(LOG_FILE, buffer));
	#endif

	char sam[MAX_BUFFER],
		 gtf[MAX_BUFFER],
		 new_sam[MAX_BUFFER],
		 mode[MAX_BUFFER];
	parse_opt(argc, argv, gtf, sam, new_sam, mode);

	if (retain_junk)
		E("-j specified; SAM file will contain junk reads!\n");

	zaman_last();


	E("Parsing GTF file %s ...\n", realpath(gtf, buffer));
	ga.parse_gtf(gtf);
	E("done in %d seconds!\n", zaman_last());

	E("Parsing SAM file %s ...\n", realpath(sam, buffer));
	parse_sam(sam);
	E("done in %d seconds!\n", zaman_last());

//	foreach (pt1, PTs_set)
//		L("%s\n", string(pt1->transcript->gene->name /*+ "." + pt1->transcript->name*/ + pt1->signature).c_str());	
//	L("======\n");

//	foreach (p, PT_single_count) {
//		PTsp pt1 = p->first.first;
///		L("%s + ", string(pt1->transcript->gene->name /*+ "." + pt1->transcript->name*/ + pt1->signature).c_str());	
	//	pt1 = p->first.second;
	//	L("%s\n", string(pt1->transcript->gene->name /*+ "." + pt1->transcript->name*/ + pt1->signature).c_str());	
//
//	}
//	L("======\n");

#ifndef LOGIFY_ONLY
	do_orman(ga, reads, read_length);

	E("Writing result to %s ...\n", new_sam);
	resolve();
	E("Resolve done..., %d sec\n", zaman_last());
	write_sam(sam, new_sam);
	E("done in %d seconds!\n", zaman_last());
#endif
#ifdef LOGIFY
	fclose(flog);
#endif

	return 0;
}

