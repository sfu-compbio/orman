/// 786

#include "common.h"
#include "rescue.h"

using namespace std;

/***********************************************************************/

void do_single (const genome_annotation &ga, vector<struct read> &reads) {
	/*E("Discarding multi-reads ...\n");
	foreach (ri, reads) {
		struct read &r = ri->second;
		if (r.entries.size() > 1) 	 // multi read
			r.entries.clear();
	}
	E("done in %d seconds!\n", zaman_last());*/
}

/***********************************************************************/

void do_rescue (const genome_annotation &ga, vector<struct read> &reads) {
	E("Initializing RESCUE ...\n");
	vector<int> transcript_unique(ga.get_transcript_count(), 0);
	// find uniqs
	foreach (ri, reads) {
		struct read &r = *ri;
		if (r.entries.size() == 1) {	 // single read
			auto *t = r.entries.begin()->partial.first->transcript;
			transcript_unique[t->id]++;
		}
	}
	E("done in %d seconds!\n", zaman_last());

	// do rescue
	E("Applying RESCUE strategy ...\n");
	foreach (ri, reads) {
		struct read &r = *ri;
		if (r.entries.size() < 2)  
			continue;

		double sum = 0;
		foreach (ti, r.entries)
			sum += transcript_unique[ti->partial.first->transcript->id] 
						/ ti->partial.first->transcript->length();
		sum = ceil(sum);

		auto it = r.entries.begin(); // result
		if (int(sum) > 0) {
			int ra = rand() % int(sum);
			sum = 0;
			foreach (ti, r.entries) {
				sum += transcript_unique[ti->partial.first->transcript->id] 
							/ ti->partial.first->transcript->length();
				if (ra < sum) {
					it = ti;
					break;
				}
			}
		}

		auto rhs = *it;
		r.entries.clear();
		r.entries.push_back(rhs);
	}
	E("done in %d seconds!\n", zaman_last());
}

