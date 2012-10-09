/// 786

#include "annotation.h"


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

template<typename T>
T &get_or_insert(map<string, T> &m, const string &k) {
	T &t = m[k];
	return t;
}

void genome_annotation::parse_gtf (const char *gtf_file) {
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
			auto it = chromosomes.find(gtf_seqname);
			int chromosome;
			if (it != chromosomes.end()) 
				chromosome = it->second;
			else { 
				chromosome = chromosomes.size();
				chromosomes.insert(make_pair(gtf_seqname, chromosome));
			}

			parse_gtf_attributes(buffer + of, attributes);

			auto &g = get_or_insert(genes, attributes["gene_id"]);
			transcript &t = get_or_insert(g.transcripts, attributes["transcript_id"]);
			t.exons.push_back(exon(gtf_start, gtf_end));
			
			g.chromosome = chromosome;
		}
	}
	transcript_count = 0;
	foreach (gi, genes) 
		foreach (ti, gi->second.transcripts) {
			transcript &t = ti->second;

			t.gene = &gi->second;
			t.name = ti->first;
			t.id   = transcript_count++;

			sort(t.exons.begin(), t.exons.end());
			for (int i = 0; i < t.exons.size(); i++) {
				t.exons[i].id = i;
				t.exons[i].transcript = &t;

				exon_tree.insert(interval(t.exons[i].start, t.exons[i].end), &t.exons[i]);
			}
		}

	fclose(fi);
}

void genome_annotation::get_exons (const interval &i, vector<exon*> &exons) {
	exon_tree.enumerate(i, exons);
}

genome_annotation::gene *genome_annotation::find_gene (uint32_t pos, int chr) {
	vector<exon*> ex;
	exon_tree.enumerate(pos, ex);
	
	foreach (e, ex) {
		if ((*e)->transcript->gene->chromosome == chr) 
			return (*e)->transcript->gene;
	}

	ex.clear();
	exon_tree.enumerate(pos - 5, ex);
	exon_tree.enumerate(pos + 5, ex);
	foreach (e, ex) {
		if ((*e)->transcript->gene->chromosome == chr)
			return (*e)->transcript->gene;
	}

	return 0;
}

