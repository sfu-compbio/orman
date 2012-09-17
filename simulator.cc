/// 786

#include "common.h"
using namespace std;

/***************************************************************************************/



/***************************************************************************************/

string genome;
map<string, uint32_t> chromosomes;

void parse_fasta (const char *fasta_file) {
	FILE *fi = fopen(fasta_file, "r");

	char buffer[MAX_BUFFER];
	genome = "";

	while (fgets(buffer, MAX_BUFFER, fi)) {
		buffer[strlen(buffer) - 1] = 0;
		if (buffer[0] == '>') {
			chromosomes[buffer + 1] = genome.size();
			continue;
		}
		genome += buffer;
	}
	E("\tRead %'u nucleotides\n", genome.size());
	E("\tRead %d chromosomes\n", chromosomes.size());

	fclose(fi);
}

/***************************************************************************************/

struct transcript {
	map<int, pair<int, int> > exons;
	string chromosome, gene, name;
	char strand;
	int expression;
	bool partial;

	transcript() :
		expression(0), partial(0), strand('+') {}
};
map<string, transcript>	  transcripts;
map<string, set<string> > genes;

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

	map<string, string> attributes;

	while (fgets(buffer, MAX_BUFFER, fi)) {
		int of; sscanf(buffer, "%s %s %s %u %u %c %c %c%n",
				gtf_seqname, gtf_source, gtf_feature,
				&gtf_start, &gtf_end, &gtf_score, &gtf_strand, &gtf_frame,
				&of);

		if (!strcmp(gtf_feature, "exon")) {
			parse_gtf_attributes(buffer + of, attributes);

			auto it = chromosomes.find(gtf_seqname);
			if (it == chromosomes.end()) {
			//	E("GTF: Unknown chromosome %s\n", gtf_seqname);
				continue;
			}
			
			transcript &t = transcripts[attributes["transcript_id"]];
			t.chromosome = it->first;
			t.strand = gtf_strand;
			t.gene = attributes["gene_id"];
			t.name = attributes["transcript_id"];
			genes[attributes["gene_id"]].insert(attributes["transcript_id"]);

			int exon_id = atoi(attributes["exon_number"].c_str());
			// if (exon_id >= t.exons.size())
			//	t.exons.resize(exon_id + 1);
			t.exons[exon_id] = make_pair(gtf_start, gtf_end);
		}
	}
	E("\tRead %'d transcripts\n", transcripts.size());

	fclose(fi);
}

/***************************************************************************************/

int total_reads = 0;
void parse_expression (const char *exp_file) {
	FILE *fi = fopen(exp_file, "r");

	char buffer[MAX_BUFFER], bf[MAX_BUFFER];
	int  exp,
		  count = 0;
	while (fgets(buffer, MAX_BUFFER, fi)) {
		sscanf(buffer, "%s %d", bf, &exp);
		auto it = transcripts.find(bf);
		if (it == transcripts.end()) {
			E("EXP: Transcript %s not found\n", bf);
			continue;
		}
		it->second.expression = exp;
		count++;
		total_reads += exp;
	}
	E("\tRead %'d expression levels\n", count);

	fclose(fi);
}

/***************************************************************************************/

int xoverlap(const pair<int, int> &a, const pair<int, int> &b) {
	int st = max(a.first, b.first);
	int ed = min(a.second, b.second);
	if (ed < st) return 0;
	return ed - st + 1;
}

bool check (transcript *tp, int skip_pos) {
	const set<string> &sp = genes[tp->gene];
	foreach (tn, sp) {
		transcript &tx = transcripts[*tn];

		if (tp->exons.size() - 1 != tx.exons.size()) // ok
			continue;

		int tpi = 1, txi = 1; // they start at 1..x
		int overlap = 0, total = 0;
		for (;;) {
			if (tpi == skip_pos) tpi++;
			if (tpi == tp->exons.size() + 1) break;
			int o = xoverlap(tp->exons[tpi], tx.exons[txi]);
			if (!o) {
				overlap = 0;
				break;
			}
			overlap += o;			
			total += tp->exons[tpi].second - tp->exons[tpi].first + 1;
			tpi++; txi++;
		}
		if (overlap / double(total) > 0.8)
			return 0;
	}
	return 1;
}

void select_partials (int percent = 10) {
	int selected = 0;
	map<string, transcript> new_transcripts;
	vector<transcript*> tp;
	foreach (t, transcripts) if (t->second.expression && t->second.exons.size() > 1) 
		tp.push_back(& t->second);

	int estimated = 0.9 * (percent / 100.0) * tp.size();
	int found = 0, treshold_c = 0, fails = 0;
	while (found < estimated && treshold_c < 10 * estimated) {
		treshold_c++;
		
		transcript *t = tp[ rand() % tp.size() ];
		if (t->partial) continue;
		int p = 1 + rand() % t->exons.size();

		if (check(t, p)) {
			t->exons.erase(p);
			t->partial = 1;
			found++;
		}
		else fails++;
	}
	E("\tSingle skips: estimated %'d, selected %'d, fails %'d\n", estimated, found, fails);

	estimated = 0.1 * (percent / 100.0) * tp.size();
	found = 0, treshold_c = 0, fails = 0;
	while (found < estimated && treshold_c < 10 * estimated) {
		treshold_c++;
		
		transcript *t = tp[ rand() % tp.size() ];
		if (t->partial) continue;
		int p = rand() % t->exons.size();

		if (check(t, p) && check(t, p + 1)) {
			transcript &tx = (new_transcripts[t->name + "_N"] = *t);
			tx.exons.erase(p + 1); 			tx.partial = 1;
			total_reads += tx.expression;
			found++;
	
			t->exons.erase(p);			t->partial = 1;
		}
		else fails++;
	}
	E("\tDouble skips: estimated %'d, selected %'d [%'d], fails %'d\n", estimated, found, new_transcripts.size(), fails);

	foreach (t, new_transcripts)
		transcripts[t->first] = t->second;
}

#define P(x) ((x) == 3 ? 'T' : ((x) == 2 ? 'G' : ((x) == 1 ? 'C' : 'A')))
int error_reads = 0;
int error_limit = 0;
int error_position[MAX_BUFFER];

inline int ohm (const string &src, char *dest, int s, int l) {
	int Px;
	for (int i = 0; i < l; i++)
		dest[i] = toupper(src[s + i]);
	dest[l] = 0;
	for (Px = 0; Px < l && (rand() % total_reads) < error_limit; Px++) {
		int r = error_position[rand() % 100];
		int n = rand() % 4;
		if (P(n) == dest[r]) n = (n + 1) % 4;
		dest[r] = P(n);
	}
	return Px;
}

int RC = 0;
void generate (const string &t_id, const transcript &t, int read_size) {
	string sequence = "";
	uint32_t offset = chromosomes[t.chromosome];

	vector<int> positions;
	foreach (e, t.exons) {
		sequence += genome.substr(e->second.first + offset, e->second.second - e->second.first + 1);
		for (int i = e->second.first; i <= e->second.second; i++)
			positions.push_back(i);
	}

	if (sequence.size() < read_size)
		return;

	int length = sequence.size() - read_size + 1;
	int jump  = max(1, length / t.expression);
	int count = max(1, t.expression / (length / jump));

	char tx[MAX_BUFFER];
	int generated = 0;
	for (int i = 0; i < length; i += jump) {
		for (int c = 0; c < count; c++) {

			int st = i + rand() % min(jump, length - i);
			int ch = ohm(sequence, tx, st, read_size);
			L(">SIM.%d chr %2s transcript %s read %d position %d length %d error %d\n", 
					RC++, t.chromosome.c_str(), t_id.c_str(), 
					generated, positions[st], read_size, ch);
			L("%s\n", tx);

			error_reads += (ch != 0);

			generated++;
			if (generated == t.expression) goto end;
		}
	}
end:
	while (generated < t.expression) {
		int st = rand() % length;
		int ch = ohm(sequence, tx, st, read_size);
		L(">SIM.%d chr %2s transcript %s read %d position %d length %d error %d\n", 
				RC++, t.chromosome.c_str(), t_id.c_str(), 
				generated, positions[st], read_size, ch);
		L("%s\n", tx);
		error_reads += (ch != 0);
		generated++;
	}
}

void generate (int read_size = 75, int novel_percentage = 10, int error_percentage = 1) {
	select_partials(novel_percentage);
	error_limit = (error_percentage / 100.0) * total_reads; 

	double curve[MAX_BUFFER];
	for (int x = 0; x < 100; x++) {
		curve[x] = - 8.34987e-14 * pow(x, 7) + 3.00017e-11 * pow(x, 6) - 4.30349e-9 * pow(x, 5) + 
			3.13868e-7 * pow(x, 4) - 0.0000121402 * pow(x, 3) + 0.000235277 * pow(x, 2) - 0.0019125 * x + 0.01;
		curve[x] *= 100.0 * 0.94715;
		if (x)
			curve[x] += curve[x - 1];
	}
	for (int i = 0; i < 100; i++) {
		int p = i;
		while (p && i < ceil(curve[p])) p--;
		while (p && i > ceil(curve[p])) p++;
		error_position[i] = p * (read_size / 100.0);
	}
	foreach(t, transcripts) if (t->second.expression) { 
		generate(t->first, t->second, read_size);
	}
	E("\tError reads: est %d generated %d\n", error_limit, error_reads);
}

/***************************************************************************************/
/* SKIP illumina model 20% sides
 *
 *
 * INDEL
 * 200 exons >= 15, introduce er 6-flank-6
	*     80%  3-indel
 */

int main (int argc, char **argv) {
	setlocale(LC_ALL, "");
	srand(time(0));
	E("Usage: simulator 1:Reference-FASTA 2:Reference-GTF 3:Reference-Expression-Levels 4:Result-GTF\nResult Output: stdout\n");

	
	char buffer[MAX_BUFFER];
	zaman_last();
	E("Parsing FASTA file %s ...\n", realpath(argv[1], buffer));
	parse_fasta(argv[1]);
	E("done in %d seconds!\n", zaman_last());

	E("Parsing GTF file %s ...\n", realpath(argv[2], buffer));
	parse_gtf(argv[2]);
	E("done in %d seconds!\n", zaman_last());

	E("Parsing expression file %s ...\n", realpath(argv[3], buffer));
	parse_expression(argv[3]);
	E("done in %d seconds!\n", zaman_last());

	E("Simulating ...\n");
	generate();
	E("done in %d seconds!\n", zaman_last());

	FILE *fo = fopen(argv[4] /* "out.gtf" */, "w");
	foreach (t, transcripts) {
		foreach (e, t->second.exons) 
			fprintf(fo, "%s protein_coding exon %d %d . %c . gene_id \"%s\"; transcript_id \"%s\"; exon_id \"%d\" novel \"%s\";\n",
					t->second.chromosome.c_str(), e->second.first, e->second.second, t->second.strand,
					t->second.gene.c_str(), t->first.c_str(), e->first, t->second.partial ? "yes" : "no");
	}
	fclose(fo);

	return 0;
}
