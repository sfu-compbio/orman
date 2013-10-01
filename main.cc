/* 786 
 *
 * Copyright (c) 2012, 2013, Simon Fraser University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the name of the Simon Fraser University nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author         : Ibrahim Numanagic
 * Email          : inumanag AT sfu DOT ca
 * Last Update    : 30. ix 2013.
 */

#include "common.h"
#include "annotation.h"
#include "interval.h"
#include "partial.h"
#include "orman.h"
#include "parser.h"
using namespace std;

#ifdef LOGIFY_ONLY
#define LOGIFY
#endif

#ifdef LOGIFY
const char *DEBUG_LOG_FILE 	= "run/orman.dbg";
const char *LOG_FILE = "run/orman.log";
FILE *flog;
#define LOG(c,...) fprintf(flog,c,##__VA_ARGS__)
#endif

int optThreads = 1;
int read_length = 0;
const char *gene_sam_flag    = "YG";
const char *partial_sam_flag = "YP";

/*******************************************************************************/

genome_annotation ga;

/*******************************************************************************/

void resolve (void) {
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
	sort(crappy_reads.begin(), crappy_reads.end());
	
	for (int i = 1; i < single_maps.size(); i++) 
		assert(single_maps[i].first != single_maps[i - 1].first);		
}

void write_sam (const char *old_sam, const char *new_sam) {
	FILE *fi = fopen(old_sam, "r"),
		 *fo = fopen(new_sam, "w");

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
		else if (j < crappy_reads.size() && line == crappy_reads[j]) {
			int l = strlen(buffer) - 1;
			sprintf(buffer + l, "\t%s:Z:_\n", gene_sam_flag);
			fputs(buffer, fo);
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
}

/*******************************************************************************/

bool file_exists (const string &s) {
	return access(s.c_str(), F_OK) != -1;
}

void parse_opt (int argc, char **argv, char *gtf, char *sam, char *newsam) {
	int opt;
	struct option long_opt[] = {
		{ "help",    0, NULL, 'h' },
		{ "gtf",  	 1, NULL, 'g' },
		{ "sam",	 1, NULL, 's' },
		{ "threads", 1, NULL, 't' },
		{ NULL,     0, NULL,  0  }
	};
	const char *short_opt = "hg:s:t:";
	do {
		opt = getopt_long (argc, argv, short_opt, long_opt, NULL);
		switch (opt) {
			case 'h':
				E("Usage: orman -t [thread count] -g [gtf file] -s [sam file] [output file]");
				E("Please check README.txt or visit http://orman.sf.net for the usage explanation.");
				exit(0);
			case 'g':
				strncpy(gtf, optarg, MAX_BUFFER);
				break;
			case 't':
				optThreads = min(optThreads, atoi(optarg));
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
	if (optind < argc)
		strncpy(newsam, argv[optind], MAX_BUFFER);
	else
		throw string("No output file specified");

	if (optThreads <= 0)
		throw string("Invalid thread count");
	if (!file_exists(gtf))
		throw string("GTF file does not exist");
	if (!file_exists(sam))
		throw string("SAM file does not exist");
}

int main (int argc, char **argv) {
	setlocale(LC_ALL, "");
	char buffer[MAX_BUFFER];

	E("ORMAN v1.1 (C) 2013 Simon Fraser University. All rights reserved.\n");
	#ifdef LOGIFY
		E("\tLog status: enabled\n");
		flog = fopen(LOG_FILE, "w");
		E("\tLog file: %s\n", realpath(LOG_FILE, buffer));
	#endif

	char sam[MAX_BUFFER] = {0};
	char gtf[MAX_BUFFER] = {0};
	char new_sam[MAX_BUFFER] = {0};

	try {
		optThreads = (int)sysconf(_SC_NPROCESSORS_ONLN) - 1;
		parse_opt(argc, argv, gtf, sam, new_sam);
		E("Using %d threads\n", optThreads);

		zaman_last();

		E("Parsing GTF file %s ...\n", realpath(gtf, buffer));
		ga.parse_gtf(gtf);
		E("done in %d seconds!\n", zaman_last());

		E("Parsing SAM file %s ...\n", realpath(sam, buffer));
		parse_sam(sam);
		E("done in %d seconds!\n", zaman_last());

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
	}
	catch (string &ex) {
		E("Error: %s\n", ex.c_str());
		exit(1);
	}
	catch (...) {
		E("Unknown error!\n");
		exit(1);
	}

	return 0;
}

