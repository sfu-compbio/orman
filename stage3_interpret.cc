/// 786

#include "common.h"
#include "interval.h"

using namespace std;


struct pt_t {
	int t; int p;
	char used;


	pt_t(): used(0) {}
	pt_t(int t, int p): t(t), p(p), used(0) {}
} ;
map<string, pt_t > result;

void read_stage2 (const char *stage2_file) {
	FILE *fi = fopen(stage2_file, "r");

	char buffer[MAX_BUFFER];
	int t, p;

	while (fscanf(fi, "%s %d", buffer, &t) != EOF) {
		if (t == -1) p = -1;
		else fscanf(fi, "%d", &p);
		result[buffer] = pt_t(t, p);
	}
}

void modify_sam (const char *sam_file, const char *stage1_file, const char *out_file) {
	FILE *fi = fopen(sam_file, "r"),
		  *f1 = fopen(stage1_file, "r");
	FILE *fo = fopen(out_file, "w");

	char buffer[MAX_BUFFER],
		  bf[MAX_BUFFER],
		  sam_name[MAX_BUFFER];

	int nl = 0, L=0;
	while (fgets(buffer, MAX_BUFFER, fi)) {
		if (buffer[0] == '@') {
			fputs(buffer, fo);
			continue;
		}
		fgets(bf, MAX_BUFFER, f1);

		char *cx;
		if ((cx = strstr(buffer, "NH:i:1")) != 0) {
			if (isdigit(cx[6]))
			{
				//E("%s\n", cx);
				//continue;
			}
			else {
				fputs(buffer, fo);
				continue;
			}
		}

		int of = 0, o, ni;
		sscanf(bf, "%s %d%n", sam_name, &ni, &o); of += o;
		
		auto &ptr = result[sam_name];
		if (ptr.used) continue;
		while (ni--) {
			int t, p; sscanf(bf + of, "%d %d%n", &t, &p, &o); of += o;
			if (t != -1 && ptr.t == t && ptr.p == p)  {
				fputs(buffer, fo);
				nl++;
				ptr.used = 1;
				break;
			}
		}
	//	E("\b\b\b\b\b\b\b\b\b\b%10d", L++);
	}
	E("\twritten %'d lines\n", nl);

	fclose(fi);
	fclose(f1);
	fclose(fo);
}

int main (int argc, char **argv) {
	setlocale(LC_ALL, "");
	E("Usage: exec SAM stage1 stage2 outSAM\n");

	char buffer[MAX_BUFFER];
	zaman_last();
	E("Parsing stage2 file %s ... ", realpath(argv[3], buffer));
	read_stage2(argv[3]);
	E("done in %d seconds!\n", zaman_last());

	E("Modifying SAM file %s\n\tstage1 %s ...\n          ", argv[1], realpath(argv[2], buffer));
	modify_sam(argv[1], argv[2], argv[4]);
	E("done in %d seconds!\n", zaman_last());


	return 0;
}

