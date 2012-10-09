/// 786

#include "common.h"
#include "annotation.h"
#include "interval.h"

using namespace std;

const char *DEBUG_LOG_FILE = "uniqorn.dbg";
const char *LOG_FILE = "uniqorn.log";

typedef genome_annotation::transcript transcript;
typedef genome_annotation::exon exon;

genome_annotation ga;

struct indel {
	int id;
	uint32_t start;
	int length;
	string insert;
	string line;

	vector<pair<int, string> > 
		support, supportN;

	indel (void) {}
	indel (int i, uint32_t s, int l, const string &in, const string &li) :
		id(i), start(s), length(l), insert(in), line(li) {}
};
vector<map<uint32_t, indel> > 
	insertions, deletions;

/*******************************************************************************/

bool parse_cigar (int chromosome, uint32_t start_pos, const char *cigar, const char *read, int line, const char *sam_line) {
	int num = 0;
	while (*cigar) {
		if (isdigit(*cigar))
			num = 10 * num + (*cigar - '0');
		else {
			if (*cigar == 'M') {
		//		result.push_back(interval(start_pos, start_pos + num));
				start_pos += num;
				read += num;
			}
			else if (*cigar == 'D') { // delete FROM REFERENCE GENOME
			//	indels.push_back(indel(start_pos, start_pos + num, indel::indel_type::DELETE));
			//	start_pos += num;
				auto i = deletions[chromosome].find(start_pos);
				if (i != deletions[chromosome].end() && i->second.length == num) 
					i->second.support.push_back(make_pair(line, sam_line));
			}
			else if (*cigar == 'N') // TODO check is the same as DEL
				start_pos += num;
			else if (*cigar == 'I') { // insert IN READ 
			//	indel x(start_pos, start_pos + num, indel::indel_type::INSERT);
			//	x.insert = string(read, num);
			//	indels.push_back(x);
				auto i = insertions[chromosome].find(start_pos);
				if (i != insertions[chromosome].end() && i->second.length == num) { 
					if (string(read, num) == i->second.insert)
						i->second.support.push_back(make_pair(line, sam_line));
					else 
						i->second.supportN.push_back(make_pair(line, sam_line));
				}

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

void parse_vcf (const char *vcf_file) {
	FILE *fi = fopen(vcf_file, "r");

	char 		buffer[MAX_BUFFER];
	char 		vcf_chromosome[MAX_BUFFER],
		  		vcf_id[MAX_BUFFER],
		  		vcf_A[MAX_BUFFER],
				vcf_B[MAX_BUFFER];
	uint32_t vcf_pos;

	int line = 0;
	while (fgets(buffer, MAX_BUFFER, fi)) {
		sscanf(buffer, "%s %u %s %s %s", vcf_chromosome, &vcf_pos, vcf_id, vcf_A, vcf_B);
		
		int chr = ga.get_chromosome(vcf_chromosome);
		if (chr == -1) {
			E("Invalid VCF line %d: %s\n", line, buffer);
			continue;
		}
		if (chr >= insertions.size())
			insertions.resize(chr + 1);
		if (chr >= deletions.size())
			deletions.resize(chr + 1);

		if (chr != -1) {
			vcf_pos++;
			int len = strlen(vcf_A) - strlen(vcf_B);

			if (len < 0) // insert
				insertions[chr][vcf_pos] = indel(line, vcf_pos, -len, vcf_A + 1, buffer);
			else 
				deletions[chr][vcf_pos] = indel(line, vcf_pos, len, "", buffer);
		}
		line++;
	}

	fclose(fi);

}

/*******************************************************************************/

string fixname (const char *c, uint32_t sam_flag) {
//	return string(c);
	return string(c) + "/" + string((sam_flag & 0x40) ? "1" : "2");
}

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

	int line = 0;
	while (fgets(buffer, MAX_BUFFER, fi)) {
		if (buffer[0] == '@') 
			continue;

		sscanf(buffer, "%s %u %s %u %u %s %s %d %d %s", 
				sam_name, &sam_flag, sam_rname, &sam_pos, &sam_mapq, sam_cigar, sam_rnext, &sam_pnext, &sam_tlen, sam_read);
		string name = fixname(sam_name, sam_flag);
		int chr = ga.get_chromosome(sam_rname);
		if (chr != -1) 
			parse_cigar(chr, sam_pos, sam_cigar, sam_read, line, buffer);
		line++;
	}

	fclose(fi);
}
/*******************************************************************************/

int main (int argc, char **argv) {
	setlocale(LC_ALL, "");
	zaman_last();

	char buffer[MAX_BUFFER];
	
	E("Parsing GTF file %s ...\n", realpath(argv[1], buffer));
	ga.parse_gtf(argv[1]);
	E("done in %d seconds!\n", zaman_last());

	E("Parsing VCF file %s ...\n", realpath(argv[2], buffer));
	parse_vcf(argv[2]);
	E("done in %d seconds!\n", zaman_last());


	E("Parsing SAM file %s ...\n", realpath(argv[3], buffer));
	parse_sam(argv[3]);
	E("done in %d seconds!\n", zaman_last());


	int ok = 0;
	foreach (I, deletions) foreach (i, *I)
		if (i->second.support.size()) {
			L("VCF Line %d support %d:\n", i->second.id, i->second.support.size());
			L("\t%s\n", i->second.line.c_str());
			foreach (s, i->second.support) L("\t%20d %s\n", s->first, s->second.c_str());
			ok++;
		}
	E("Expressed deletions %d / %d\n", ok, deletions.size());

	ok = 0;
	foreach (I, insertions) foreach (i, *I)
		if (i->second.support.size()) {
			L("VCF Line %d support %d:\n", i->second.id, i->second.support.size());
			L("\t%s\n", i->second.line.c_str());
			foreach (s, i->second.support) L("\t%20d %s\n", s->first, s->second.c_str());
			ok++;
		}
	E("Expressed insertions %d / %d\n", ok, insertions.size());

	ok = 0;
	foreach (I, insertions) foreach (i, *I)
		if (!i->second.support.size() && i->second.supportN.size()) {
			L("VCF Line %d supportN %d:\n", i->second.id, i->second.supportN.size());
			L("\t%s\n", i->second.line.c_str());
			foreach (s, i->second.supportN) L("\t%20d %s\n", s->first, s->second.c_str());
			ok++;
		}
	E("Almost expressed insertions %d / %d\n", ok, insertions.size());

	return 0;
}

