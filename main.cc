/// 786

#include "common.h"
#include "annotation.h"
#include "interval.h"
#include "partial.h"
#include "orman.h"
#include "parser.h"
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
const char *gene_sam_flag    = "YG";
const char *partial_sam_flag = "YP";

FILE *flog;
#define LOG(c,...) fprintf(flog,c,##__VA_ARGS__)

/*******************************************************************************/

genome_annotation ga;

/*******************************************************************************/

void resolve () {
	foreach (r, reads) {
		if (r->entries.size() == 1) {
			assert(r->entries.begin()->line.first != -1);
			single_maps.push_back(make_pair(
						r->entries.begin()->line.first, 
						r->entries.begin()->partial->get_gene()));
			if (r->entries.begin()->line.second != -1)
				single_maps.push_back(make_pair(
						r->entries.begin()->line.second, 
						r->entries.begin()->partial->get_gene()));
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

//	foreach (p, PT_single_count) {
//		PTsp pt1 = p->first.first;
///		L("%s + ", string(pt1->transcript->gene->name /*+ "." + pt1->transcript->name*/ + pt1->signature).c_str());	
	//	pt1 = p->first.second;
	//	L("%s\n", string(pt1->transcript->gene->name /*+ "." + pt1->transcript->name*/ + pt1->signature).c_str());	
//	}
//	L("======\n");

#ifndef LOGIFY_ONLY
	do_orman(ga, reads, read_length);

	E("Writing result to %s ...\n", new_sam);
	resolve();
	write_sam(sam, new_sam);
	E("done in %d seconds!\n", zaman_last());
#endif
#ifdef LOGIFY
	fclose(flog);
#endif

	return 0;
}

