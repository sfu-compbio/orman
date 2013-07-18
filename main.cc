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

const char *DEBUG_LOG_FILE = "orman.dbg";
const char *LOG_FILE 		= "orman.log";

char junk_file_path[MAX_BUFFER];

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

bool 		retain_junk 	= false;
int 		read_length 	= 0;
const char*		gene_sam_flag 	= "YG";

/*******************************************************************************/

set<PTs> 							PTs_set;
map<pair<PTsp, PTsp>, int> 	PT_single_count;
vector<struct read> 				reads;
genome_annotation 				ga;

int get_single_coverage(const PT &p) {
	auto i = PT_single_count.find(p);
	if (i == PT_single_count.end())
		return 0;
	else 
		return i->second;
}

/*******************************************************************************/

bool parse_cigar (uint32_t start_pos, const char *cigar, const char *read, vector<interval> &result, vector<indel> &indels) {
	int num = 0;
	
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

	if (exons.size() == 0 || parts.size() == 0) {
#ifdef LOGIFY
		get_partial_transcript_error += S("exon_span=%d; cigar_parts=%d; ", exons.size(), parts.size());
#endif
		return 0;
	}
	// disregard reads with mapping chunk <4
	foreach (pi, parts)
		if (pi->second - pi->first + 1 < 4) {
#ifdef LOGIFY
			get_partial_transcript_error += S("part_sz<4; ");
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
			get_partial_transcript_error += S("intron_beginning_exon_%d; ", e->id);
#endif
			return 0; 
		}
		// find internal skips
		int prev = parts[pi].first;
		while (pi < parts.size() && ((parts[pi].second <= e->end) || (parts[pi].second - e->end < 4))) {
			// internal exon skip
			if (prev < parts[pi].first) {
#ifdef LOGIFY
			get_partial_transcript_error += S("inexon_skip_%d; ", e->id);
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
			get_partial_transcript_error += S("intron_ending_exon_%d; ", e->id);
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

void makecand (const vector<interval> &parts, int chromosome, map<transcript*, vector<exon*> > &candidates) {
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

/*
 * get all pairs of partial transcripts for each paired read
 */
void parse_read (int64_t line1, uint32_t start_pos1, int chromosome1, const char *cigar1, const char *read1,
					  int64_t line2, uint32_t start_pos2, int chromosome2, const char *cigar2, const char *read2,
					  set<pair<tuple<PTsp, int, int64_t>, tuple<PTsp, int, int64_t>>> &result) {
	vector<interval> parts1,  parts2;
	vector<indel>    indels1, indels2;

	if (line1 >= 0) parse_cigar(start_pos1, cigar1, read1, parts1, indels1);
	if (line2 >= 0) parse_cigar(start_pos2, cigar2, read2, parts2, indels2);

	map<transcript*, vector<exon*> > candidates1, candidates2;
	trimx(parts1); makecand(parts1, chromosome1, candidates1);
	trimx(parts2); makecand(parts2, chromosome2, candidates2);

#ifdef LOGIFY
	LOG("%d & %d\n", line1, line2);
	if (candidates1.size() == 0) 
		LOG("\t/1.%09d %s no_valid_candidiates;\n", line1, cigar1);
	if (line2 != 0 && candidates2.size() == 0) 
		LOG("\t/2.%09d %s no_valid_candidiates;\n", line2, cigar2);
#endif

#ifdef LOGIFY_ONLY
	foreach (ci1, candidates1) {
		int starting_position1 = 0; 
		PTsp pt1 = get_partial_transcript(ci1->first, ci1->second, parts1, indels1, starting_position1);
		LOG("\t/1.%09d %s %s\n", line1, cigar1, pt1 
				? string(pt1->transcript->gene->name + "_" +pt1->transcript->name + "_" + pt1->signature).c_str()
				: get_partial_transcript_error.c_str());
	}
	foreach (ci2, candidates2) {
		int starting_position2 = 0; 
		PTsp pt2 = get_partial_transcript(ci2->first, ci2->second, parts2, indels2, starting_position2);
		LOG("\t/2.%09d %s %s\n", line2, cigar2, pt2 
				? string(pt2->transcript->gene->name + "_" +pt2->transcript->name + "_" + pt2->signature).c_str()
				: get_partial_transcript_error.c_str());
	}
#else
	foreach (ci1, candidates1) {
		int starting_position1 = 0; 
		PTsp pt1 = get_partial_transcript(ci1->first, ci1->second, parts1, indels1, starting_position1);
#ifdef LOGIFY
		LOG("\t/1.%09d %s %s\n", line1, cigar1, pt1 
				? string(pt1->transcript->gene->name + "_" +pt1->transcript->name + "_" + pt1->signature).c_str()
				: get_partial_transcript_error.c_str());
#endif

		// first mate HAS TO BE VALID!
		if (pt1) {
			// special case if there is no valid pair
			bool with_null = false;
			if (candidates2.size()) foreach (ci2, candidates2) {
				int starting_position2 = 0;
				PTsp pt2 = get_partial_transcript(ci2->first, ci2->second, parts2, indels2, starting_position2);
				if (pt2) 
					result.insert(make_pair(make_tuple(pt1, starting_position1, line1), make_tuple(pt2, starting_position2, line2)));
#ifdef LOGIFY
				LOG("\t/2.%09d %s %s\n", line2, cigar2, pt2 
						? string(pt2->transcript->gene->name + "_" +pt2->transcript->name + "_" + pt2->signature).c_str()
						: get_partial_transcript_error.c_str());
#endif

				with_null |= (pt2 != 0);
			}
			if (with_null)
				result.insert(make_pair(make_tuple(pt1, starting_position1, line1), make_tuple((PTsp)0, 0, line2)));
		}
	}
#endif
}


/*******************************************************************************/

vector<pair<int64_t, genome_annotation::gene*>> single_maps;
vector<int64_t> crappy_reads;

void parse_sam (const char *sam_file) {
	FILE *fi = fopen(sam_file, "r");

	// obtain the file size
	fseek(fi, 0, SEEK_END);
	int64_t f_size = ftell(fi);
	fseek(fi, 0, SEEK_SET);

	// define storage for the same reads
	// key:   <<chr, pos>, <nextchr, nextpos>, frag_len, is_secondary_mapping>
	typedef tuple<	pair<int, uint32_t>,
			  			pair<int, uint32_t>,
						int, bool > idx_key;
	// value: <<line1, read1, cigar1>, <line2, read2 cigar2>
	typedef pair< tuple<int64_t, string, string>, 
			  		  tuple<int64_t, string, string> > idx_val;
	// UGLY! but works for now ...
	multimap<idx_key, idx_val> idx;

	char 		buffer[MAX_BUFFER];
	char 		sam_name[MAX_BUFFER],
		  		sam_rname[MAX_BUFFER],
		  		sam_cigar[MAX_BUFFER],
				sam_read[MAX_BUFFER],
				sam_rnext[MAX_BUFFER];
	uint32_t sam_flag, sam_pos;
	int32_t  sam_pnext, sam_tlen;
	uint8_t  sam_mapq;

	int64_t line = 0;
	int read_id = 0;
	string prev_name = "";

	E("\t%5s %15s %15s %15s\n", "%%", "Partials", "Reads", "SAM lines");
	while (fgets(buffer, MAX_BUFFER, fi)) {
		if (buffer[0] == '@') 
			continue;

		sscanf(buffer, "%s %u %s %u %u %s %s %d %d %s", 
				sam_name, &sam_flag, sam_rname, &sam_pos, &sam_mapq, sam_cigar, sam_rnext, &sam_pnext, &sam_tlen, sam_read);

		if (prev_name != string(sam_name) || feof(fi))
		{
			// iterate through all reads and obtain all pairs of partial transcripts
			set<pair<tuple<PTsp, int, int64_t>, tuple<PTsp, int, int64_t>>> result, tmp;
#ifdef LOGIFY
			LOG("%s =%d\n", prev_name.c_str(), idx.size());
#endif
			foreach (i, idx) {
				auto &r1 = i->second.first,
					  &r2 = i->second.second;
				parse_read(
					get<0>(r1), get<0>(i->first).second, get<0>(i->first).first, get<2>(r1).c_str(), get<1>(r1).c_str(),
					get<0>(r2), get<1>(i->first).second, get<1>(i->first).first, get<2>(r2).c_str(), get<1>(r2).c_str(),
					tmp
				);
				if (retain_junk && tmp.size() == 0) { // crappy?!
					if (get<0>(r1) >= 0) crappy_reads.push_back(get<0>(r1));
					if (get<0>(r2) >= 0) crappy_reads.push_back(get<0>(r2));
				}
				result.insert(tmp.begin(), tmp.end());
			}

#ifndef LOGIFY_ONLY
			// process ONLY if we have multi-mappings!
			if (result.size() > 1) { 
				if (read_id >= reads.size())
					reads.resize(read_id + 10000);
				foreach (i, result) {
					reads[read_id].entries.push_back(read::read_entry(
						make_pair(get<0>(i->first), get<0>(i->second)), // PT 
						make_pair(get<2>(i->first), get<2>(i->second)), // line 
						make_pair(get<1>(i->first), get<1>(i->second))  // PT start pos
					));
				}
				read_id++;
			}
			// otherwise, just update the single-mapping partial counter
			else if (result.size() == 1) {
				auto k = make_pair(get<0>(result.begin()->first), get<0>(result.begin()->second));
				auto i = PT_single_count.find(k);
				if (i != PT_single_count.end()) 
					i->second++;
				else 
					PT_single_count[k] = 1;
				single_maps.push_back(make_pair(get<2>(result.begin()->first),  k.first->get_gene()));
				single_maps.push_back(make_pair(get<2>(result.begin()->second), k.second->get_gene()));
			}
			// if there are no valid PTs, discard
			else ;
#endif

			prev_name = string(sam_name);
			idx.clear();
		}
		// fix chromosome value
		if (string(sam_rnext) == "=")
			strcpy(sam_rnext, sam_rname);
		// get chromosome values
		int chr1 = ga.get_chromosome(sam_rname),
			 chr2 = ga.get_chromosome(sam_rnext);
		// check chromosomes
		if (chr1 == -1 || chr2 == -1) { 
			// do not process
			line++; continue; 
		}

		// is it first mate?
		if ((sam_flag & 0x8) || (sam_flag & 0x40)) {
			auto k = idx_key(make_pair(chr1, sam_pos), 
					make_pair(chr2, sam_pnext),
					abs(sam_tlen), sam_flag & 0x100); // not primary alignment
			auto i = idx.find(k);

			// value: <<line1, line2>, <read1, read2>, <cigar1, cigar2>
			if (i != idx.end()) {
				get<0>(i->second) = make_tuple(line, string(sam_read), string(sam_cigar));					
			}
			else
				idx.insert(make_pair(k, idx_val( make_tuple(line, string(sam_read), string(sam_cigar)), make_tuple(int64_t(-1), string(""), string("")) )));	
		}
		else {
			auto k = idx_key(make_pair(chr2, sam_pnext), 
					make_pair(chr1, sam_pos),
					abs(sam_tlen), sam_flag & 0x100);
			auto i = idx.find(k);

			// value: <<line1, line2>, <read1, read2>, <cigar1, cigar2>
			if (i != idx.end()) {
				get<1>(i->second) = make_tuple(line, string(sam_read), string(sam_cigar));					
			}
			else
				idx.insert(make_pair(k, idx_val( make_tuple(int64_t(-1), string(""), string("")), make_tuple(line, string(sam_read), string(sam_cigar)) )));	
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
			E("OOOOOOOOOOOOOps!");
			exit(1);
		}
	}

	fclose(fi);
}

void resolve () {
	foreach (r, reads) {
		if (r->entries.size() == 1) {
				single_maps.push_back(make_pair(
							r->entries.begin()->line.first, 
							r->entries.begin()->partial.first->get_gene()));
				if (r->entries.begin()->partial.second) single_maps.push_back(make_pair(
							r->entries.begin()->line.second, 
							r->entries.begin()->partial.second->get_gene()));
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

	E("\t%5s %15s\n", "%%", "SAM lines");
	while (fgets(buffer, MAX_BUFFER, fi)) {
		if (buffer[0] == '@') {	
			fputs(buffer, fo);
			continue;
		}
		// TODO fix NH:i:...
		if (i < single_maps.size() && line == single_maps[i].first) {
			int l = strlen(buffer) - 1;
			sprintf(buffer + l, " %s:A:%s\n", gene_sam_flag, 
					single_maps[i].second ? single_maps[i].second->name.c_str() : "");
			fputs(buffer, fo);
			i++;
			wr++;
		}
		else if (retain_junk && j < crappy_reads.size() && line == crappy_reads[j]) {
			fputs(buffer, fj);
			j++;
		}
		else di++;
		line++;
		if (line % (1<<14) == 0) 
			E("\r\t%5.2lf %'15d", 100.0 * double(ftell(fi)) / f_size, line);
	}
	E("\nOK, written %'d, discarded %'d, total %'d\n", wr, di, line);
	E("      crappy  %'d, total %'d\n", j, wr + j);

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
		{ "mode",	 1, NULL, 'm' },
		{ NULL,     0, NULL,  0  }
	};
	const char *short_opt = "hj:g:s:m:";
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
			case 'm':
				strncpy(mode, optarg, MAX_BUFFER);
				break;
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

	E("Behold! Former Uniqorn Republic of Orman is starting!\n\tUsage: orman -g[gtf] -s[sam] -m[mode:(single|rescue|orman)] [new_sam]\n");
	E("\tCompile time: %s\n", COMPILE_TIME);
#ifdef LOGIFY
	E("\tLog status: enabled\n");
	flog = fopen(LOG_FILE, "w");
#endif
	//freopen(DEBUG_LOG_FILE, "w", stdout);

	char sam[MAX_BUFFER],
		  gtf[MAX_BUFFER],
		  new_sam[MAX_BUFFER],
		  mode[MAX_BUFFER],
		  buffer[MAX_BUFFER];
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

//	foreach(i, PTs_set) {
//		L("%25s %20s %15d\n", i->transcript->name.c_str(), i->signature.c_str(), i->weight);
//	}
//	exit(0);

#ifndef LOGIFY_ONLY
	if (!strcmp(mode, "rescue"))
		do_rescue(ga, reads);
	else if (!strcmp(mode, "single"))
		do_single(ga, reads);
	else if (!strcmp(mode, "orman"))
		do_orman(ga, reads, read_length);
	else {
		E("Unknown mode %s!\n");
		exit(1);
	}

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

