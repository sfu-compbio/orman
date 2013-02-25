/// 786

#include "common.h"
#include "annotation.h"
#include "interval.h"
#include "rescue.h"
#include "partial.h"
#include "orman.h"

using namespace std;

const char *DEBUG_LOG_FILE = "orman.dbg";
const char *LOG_FILE 		= "orman.log";

bool retain_junk = false;

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

int read_length = 0;

/*******************************************************************************/

set<partial_transcript_single> partial_transcripts_single;
set<partial_transcript> partial_transcripts;
vector<struct read> reads;
genome_annotation ga;

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

const partial_transcript_single *get_partial_transcript (transcript *t, const vector<exon*> &exons, const vector<interval> &parts, const vector<indel> &indels, int &starting_position) {
	// disregard junk
	if (exons.size() == 0 || parts.size() == 0) {
		L("exsz %d pasz %d", exons.size(), parts.size());
		return 0;
	}
	// disregard reads with mapping chunk <4
	foreach (pi, parts)
		if (pi->second - pi->first + 1 < 4) {
			L("map chunk <4");
			return 0;
		}
	// disregard reads with introns and internal exon skips
	int pi = 0;
	foreach (ei, exons) {
		exon *e = *ei;
		// disregard introns
		// only allow +- 4 introns at the beginning
		if (parts[pi].first < e->start && !(e->start - parts[pi].first < 4)) { 
			L("intron1");
			return 0; 
		}
		// find internal skips
		int prev = parts[pi].first;
		while (pi < parts.size() && ((parts[pi].second <= e->end) || (parts[pi].second - e->end < 4))) {
			// internal exon skip
			if (prev < parts[pi].first) {
				L("skip");
				return 0;
			}
			prev = parts[pi].second;
			pi++;
		}
		// again, no introns
		// allow +-4 at the end as well
		if (pi < parts.size () && parts[pi].first <= e->end && !(parts[pi].second - e->end < 4)) {
			L("intron2");
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
		signature += "e" + numtostr(e->id);
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

	partial_transcript_single pt(t, signature, length, weight);
	auto pti = partial_transcripts_single.insert(pt);
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

void parse_read (int name_id, int line1, uint32_t start_pos1, int chromosome1, const char *cigar1, const char *read1,
					  					int line2, uint32_t start_pos2, int chromosome2, const char *cigar2, const char *read2) {
	vector<interval> parts1,  parts2;
	vector<indel>    indels1, indels2;
	parse_cigar(start_pos1, cigar1, read1, parts1, indels1);
	if (line2 >= 0) parse_cigar(start_pos2, cigar2, read2, parts2, indels2);

	map<transcript*, vector<exon*> > candidates1, candidates2;
	trimx(parts1); makecand(parts1, chromosome1, candidates1);
	trimx(parts2); makecand(parts2, chromosome2, candidates2);
	
	foreach (ci1, candidates1) foreach (ci2, candidates2) {
		int starting_position1 = 0, 
			 starting_position2 = 0;

	//	L("%40s ", name.c_str());
		const partial_transcript_single *pt1 = get_partial_transcript(ci1->first, ci1->second, parts1, indels1, starting_position1);
		const partial_transcript_single *pt2 = line2 >= 0 ? get_partial_transcript(ci2->first, ci2->second, parts2, indels2, starting_position2) : 0;
		
/*		if (pt) L(" %20s %20s %10d ", pt->transcript->name.c_str(), pt->signature.c_str(), starting_position); 
		L(" cigar %s tr %s ex ", cigar, ci->first->name.c_str()); foreach(e, ci->second) L("(%d,%d) ",(*e)->start, (*e)->end); L("parts "); foreach(p,parts)L("(%d %d) ", p->first, p->second); 
		L("\n");*/

	// right now strand/1 HAS TO BE VALID!
		if (pt1) { 		
		//	if (reads[name].entries.find(line)!=reads[name].entries.end())
		//		E("OOOOO! %s %d\n",name.c_str(), line);
			auto pti = partial_transcripts.insert(partial_transcript(pt1, pt2));
			const partial_transcript *pt = &(* pti.first);

			if (reads.size() <= name_id)
				reads.resize(name_id + 1000);
			reads[name_id].entries[read::read_key(line1, line2, pt)] = make_pair(
				read::read_entry(start_pos1, starting_position1), 
				read::read_entry(start_pos2, starting_position2)
			);
		}
	}
}


/*******************************************************************************/

string fixname (const char *c, uint32_t sam_flag) {
	return string(c) + "/" + string((sam_flag & 0x40) ? "1" : "2");
}

#include <tuple>
typedef pair<string, uint32_t> keyx___;
typedef std::tuple<keyx___, keyx___, int32_t, bool> key___;
struct read___ {
	int id;
	string cigar, read;
	int chr, line;
	uint32_t sam_pos, sam_flag;

	read___ () {}
	read___ (int n, int l, int sp, uint32_t sf, int c, const string &cc, const string &r) :
		id(n), line(l), sam_pos(sp), sam_flag(sf), chr(c), cigar(cc), read(r) {}
};

void parse_sam (const char *sam_file) {
	FILE *fi = fopen(sam_file, "r");

	fseek(fi, 0, SEEK_END);   // non-portable
	int64_t f_size = ftell(fi);
	fseek(fi, 0, SEEK_SET);

	// what if sam file hates us?
	multimap<key___, read___> short_index;

	char 		buffer[MAX_BUFFER];
	char 		sam_name[MAX_BUFFER],
		  		sam_rname[MAX_BUFFER],
		  		sam_cigar[MAX_BUFFER],
				sam_read[MAX_BUFFER],
				sam_rnext[MAX_BUFFER];
	uint32_t sam_flag, sam_pos;
	int32_t  sam_pnext, sam_tlen;
	uint8_t  sam_mapq;

	int line = 0;
	int read_id = -1;
	string prev_name = "";
	E("\tSAM %%\t\tPartials\t\tReads\n");
	while (fgets(buffer, MAX_BUFFER, fi)) {
		if (buffer[0] == '@') 
			continue;

		sscanf(buffer, "%s %u %s %u %u %s %s %d %d %s", 
				sam_name, &sam_flag, sam_rname, &sam_pos, &sam_mapq, sam_cigar, sam_rnext, &sam_pnext, &sam_tlen, sam_read);

		if (prev_name != string(sam_name)) {
			if (short_index.size() != 0) {
				E("File not SORTED! Line %d, read %d\n", line-1, prev_name.c_str());
				abort();
			}

			prev_name = string(sam_name);
			read_id++;
		}

		read_length = max(read_length, (int)strlen(sam_read));
		int chr = ga.get_chromosome(sam_rname);
		if (chr == -1) { line++; continue; }
		bool secondary = sam_flag & 0x100;

		if (sam_flag & 0x8) { // nema paired enda!
			parse_read(read_id, line, sam_pos, chr, sam_cigar, sam_read,
								       -1,       0,   0, sam_cigar, sam_read);
		}
		else { // uh oh
			if (string(sam_rnext) == "=")
				strcpy(sam_rnext, sam_rname);
			read___ rx(read_id, line, sam_pos, sam_flag, chr, sam_cigar, sam_read);
			auto ix = short_index.find(key___( 
						keyx___(sam_rname, sam_pos), keyx___(sam_rnext, sam_pnext), 
						abs(sam_tlen), secondary));
			
			if (ix != short_index.end()) {
				read___ r1 = ix->second, r2 = rx;
				if (r2.sam_flag & 0x40) swap(r1, r2);
				parse_read(read_id, r1.line, r1.sam_pos, r1.chr, r1.cigar.c_str(), r1.read.c_str(),
										  r2.line, r2.sam_pos, r2.chr, r2.cigar.c_str(), r2.read.c_str());
				short_index.erase(ix);
			}
			else short_index.insert(make_pair(
				key___(keyx___(sam_rnext, sam_pnext), 
					keyx___(sam_rname, sam_pos), abs(sam_tlen), secondary),
				rx));
		}

		line++;

		if (line % 8092 == 0) 
			E("\r\t%.2lf%%\t\t%'d\t\t%'d", 100.0 * double(ftell(fi)) / f_size, partial_transcripts_single.size(), reads.size());
	}
	assert(short_index.size() == 0);
	E("\n");

	fclose(fi);
}

void write_sam (const char *old_sam, const char *new_sam) {
	FILE *fi = fopen(old_sam, "r"),
		  *fo = fopen(new_sam, "w");
	FILE *fl = fopen(LOG_FILE, "w");

	char 		buffer[MAX_BUFFER];
	char 		sam_name[MAX_BUFFER],
		  		sam_rname[MAX_BUFFER],
		  		sam_cigar[MAX_BUFFER],
				sam_read[MAX_BUFFER],
				sam_rnext[MAX_BUFFER];
	uint32_t sam_flag, sam_pos;
	int32_t  sam_pnext, sam_tlen;
	uint8_t  sam_mapq;

	int wr = 0, 
		 di = 0,
		 ix = 0;
	int line = 0;

	int read_id = -1;
	string prev_name = "";

	while (fgets(buffer, MAX_BUFFER, fi)) {
		if (buffer[0] == '@') {	
			fputs(buffer, fo);
			continue;
		}
		//sscanf(buffer, "%s %u %s %u", sam_name, &sam_flag, sam_rname, &sam_pos);

		sscanf(buffer, "%s %u %s %u %u %s %s %d %d %s", 
				sam_name, &sam_flag, sam_rname, &sam_pos, &sam_mapq, sam_cigar, sam_rnext, &sam_pnext, &sam_tlen, sam_read);
		string name = sam_name; //fixname(sam_name, sam_flag);

		if (prev_name != string(sam_name)) {
			prev_name = string(sam_name);
			read_id++;
		}

		int chr = ga.get_chromosome(sam_rname);
		if (chr == -1) { 
			fprintf(fl, "%40s discarded\n", name.c_str());
			fputs(buffer, fo);
			ix++;
			line++;
			continue; 
		}


		if (read_id >= reads.size()) {
			E("Error: r_id > reads!\n");
			abort();
		}

		const auto &ri = reads[read_id];
		if (ri.entries.size() == 1 &&
			((ri.entries.begin()->first.line == line) || 
			(ri.entries.begin()->first.line2 == line)))
		{
//			const read::read_entry &re = ri->second.entries.begin()->second;

//			assert(re.chromosome == chr);
//			assert(re.position == sam_pos);

		/*	fprintf(fl, "%40s %3d %16u %20s %10d %s\n",
					name.c_str(), 
					re.chromosome, re.position,
					re.partial->transcript->name.c_str(), 
					re.partial_start,
					re.partial->signature.c_str());*/
		
			fputs(buffer, fo);
			wr++;
		}
		else if (/*ri == reads.end() ||*/ ri.entries.size() == 0) {
	//		fprintf(fl, "%40s discarded\n", name.c_str());
			if(retain_junk)
				fputs(buffer, fo);
			ix++;
		}
		else 
			di++;
		line++;
	}
	E("\tOK, written %'d, discarded %'d, forcefully retained %'d, total %'d\n", wr, di, ix, line);

	fclose(fo);
	fclose(fi);

	fclose(fl);
}

/*******************************************************************************/

void parse_opt (int argc, char **argv, char *gtf, char *sam, char *newsam, char *mode) {
	int opt; 
	struct option long_opt[] = {
		{ "help",    0, NULL, 'h' },
		{ "gtf",  	 1, NULL, 'g' },
		{ "junk",  	 1, NULL, 'j' },
		{ "sam",		 1, NULL, 's' },
		{ "mode",	 1, NULL, 'm' },
		{ NULL,     0, NULL,  0  }
	};
	const char *short_opt = "hjg:s:m:";
	do {
		opt = getopt_long (argc, argv, short_opt, long_opt, NULL);
		switch (opt) {
			case 'h':
				exit(0);
			case 'g':
				strncpy(gtf, optarg, MAX_BUFFER);
				break;
			case 'j':
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
	freopen(DEBUG_LOG_FILE, "w", stdout);

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

	if (!strcmp(mode, "rescue"))
		do_rescue(ga, reads);
	else if (!strcmp(mode, "single"))
		do_single(ga, reads);
	else if (!strcmp(mode, "orman"))
		do_orman(ga, partial_transcripts, reads, read_length);
	else {
		E("Unknown mode %s!\n");
		exit(1);
	}

	E("Writing result to %s ...\n", new_sam);
	write_sam(sam, new_sam);
	E("done in %d seconds!\n", zaman_last());
	

	return 0;
}

