/// 786

#include "common.h"
#include "orman.h"
#include <boost/heap/fibonacci_heap.hpp>
// #include <ilcplex/ilocplex.h>

using namespace std;
using namespace boost::heap;

struct cluster {
	int id;
	const partial_transcript *partial;
	// position -> [<fatreads ...>]
	map<int, set<int> > fatreads;
	bool covered;

	cluster (void) {}
	cluster (int i, const partial_transcript *p) :
		id(i), partial(p), covered(0) {}
};
vector<cluster> clusters;

struct fatread {
	// IDs of the reads
	vector<struct read*> reads;
	// Cluster ID -> position
	map<int, int> clusters;
	// Cluster ID -> HowMuch
	map<int, int> solution;
};
vector<fatread> fatreads;

/*****************************************************************************************/

void initialize_structures (const set<partial_transcript> &partials, map<string, struct read> &reads) {
	// initialize our clusters
	clusters.resize(partials.size());

	// assign clusters to partials
	int id = 0;
	map<const partial_transcript*, int> cluster_index;
	foreach (ci, partials) {
		clusters[id] = cluster(id, &*ci);
		cluster_index[&*ci] = id;
	//	E("%s.%s %d\n", ci->transcript->name.c_str(), ci->signature.c_str(), id);
		id++;
	}

	// initialize fatreads
	char buffer[MAX_BUFFER];
	map<string, int> fatread_index;
	foreach (ri, reads) {
		struct read &r = ri->second;
		if (r.entries.size() < 1)
			continue;

	//	E("%s -> ", ri->first.c_str());

		// prepare and sort signature!
		map<const partial_transcript*, int> sig;
		foreach (ei, r.entries) 
			sig.insert(make_pair(ei->second.partial, ei->second.partial_start));
		// form signature
		string signature = "";
		foreach (si, sig) {
			sprintf(buffer, "%llu+%d;", cluster_index[si->first], si->second);
	//		E("%s.%s %d ", si->first->transcript->name.c_str(), si->first->signature.c_str(), cluster_index[si->first]);
			signature += buffer;
		}
	//	E("\n");
	
		// add fatread to the index if it is not there
		auto fi = fatread_index.find(signature);
		int fid;
		if (fi == fatread_index.end()) {
			fid = fatread_index.size();
			fatread_index[signature] = fid;
			if (fatreads.size() <= fid) 
				fatreads.resize(fid + 1);
			// initialize fatread
			foreach (si, sig) {
				int cluster = cluster_index[si->first];
//				E("%d\n", cluster);
				fatreads[fid].clusters[cluster] = si->second;
				fatreads[fid].solution[cluster] = 0; // do not assign anything at the beginning
			}
		}
		else fid = fi->second;
		// assign this read to the fatread
		fatreads[fid].reads.push_back(&r);
	}
	fatread_index.clear();
	cluster_index.clear();

	// initialize clusters
	for (int fi = 0; fi < fatreads.size(); fi++) {
		assert(fatreads[fi].reads.size() > 0);
		foreach (ci, fatreads[fi].clusters) {
	//		if(ci->first==27)E("xaa\n");
			clusters[ci->first].fatreads[ci->second].insert(fi);
		}
	}

	/*int lc=0;
	foreach (c, clusters) {
		if (c->fatreads.size() == 0)
			lc++;
	}
	E(">>%d\n",lc);*/

	E("\t%'d fatreads, %'d clusters\n", fatreads.size(), clusters.size());
}

/*****************************************************************************************/

struct heap_t {
	cluster *t;

public:
	int cardinality;
	double weight;

public:
	heap_t (void): 
		t(0) {}
	heap_t (cluster *t, int c, double w): 
		t(t), cardinality(c), weight(w) {}

public:
	double value (void) const {
		return cardinality / weight;
	}
	bool operator< (const heap_t &b) const {
		return value() < b.value();
	}
};

int set_cover (int read_length) {
	fibonacci_heap<heap_t> heap;
	vector<fibonacci_heap<heap_t>::handle_type> heap_handles(clusters.size());
	for (int i = 0; i < clusters.size(); i++) {
		int covers = 0;
		// count number of reads this set covers
		foreach (pos, clusters[i].fatreads) 
			foreach (fr, pos->second) // in fatreads
				covers += fatreads[*fr].reads.size();

		// calculate number of non-covered positions!
		// and update the weight
		double weight = clusters[i].partial->weight;
		//E("[%10.2lf %5d] ", weight, clusters[i].partial->length);
		int notcovered = 0;
		int prev = 0;
		foreach (fr, clusters[i].fatreads) { // pos -> <fatreads>
			if (fr->first > prev)
				notcovered += fr->first - prev;
			prev = fr->first + read_length;
		}
		if (prev < clusters[i].partial->length)
			notcovered += clusters[i].partial->length - prev;
		if (notcovered > 1) {
		//	LE("Set %s len %d notc %d\n", t->signature.c_str(), t->length, notcovered);
		//	TODO talk with Phuong
		//	E(" (%5d / %5d) ", notcovered, clusters[i].partial->length);
			weight *= (1 + 
					10.0 * double(notcovered) / clusters[i].partial->length);
		}
	
		// E("> %15s %10s %5d %.2lf\n", clusters[i].partial->transcript->name.c_str(), clusters[i].partial->signature.c_str(), covers, weight);
		heap_handles[i] = heap.push(heap_t(&clusters[i], covers, weight));
		assert(covers > 0);
	}

	// do the greedy algorithm
	vector<char> read_visited(fatreads.size(), 0);
	int mincover = 0;
	while (!heap.empty()) {
		heap_t h = heap.top(); heap.pop();

		if (h.cardinality == 0)
			break;
		mincover++;

		foreach (pos, h.t->fatreads) 
			foreach (r, pos->second) /* in fatreads */ if (!read_visited[*r]) {
			//	fatreads[*r].solution[ h.t->id ] += fatreads[*r].reads.size(); // assign ALL to this set

				foreach (s, fatreads[*r].clusters) if (s->first != h.t->id) {
					(*heap_handles[s->first]).cardinality -= fatreads[*r].reads.size();
					heap.decrease(heap_handles[s->first]);
				}
				read_visited[*r] = 1;
			}
		h.t->covered = true;
	}



	int not_covered = 0, ncr = 0;

/*	foreach (c, clusters)
		if (!c->covered)
			not_covered++;
	E("\tNot covered %d sets\n", not_covered);
	not_covered=0;*/

	foreach (fr, fatreads) {
		char cv = 0;
		foreach (t, fr->clusters)
			cv |= clusters[t->first].covered;
		if (!cv) {
			not_covered++;
			ncr += fr->reads.size();
			//E("<<%d:\n",fr->clusters.size());
			//foreach(cx, fr->clusters)
			//	E("\t%d %d %.2lf\n", clusters[cx->first].covered, (*heap_handles[cx->first]).cardinality, (*heap_handles[cx->first]).weight);
		}
	}
	E("\t%'d non-covered fatreads (%'d reads)\n", not_covered, ncr);

	return mincover;
}

void probabilistic_assign (void) {
	foreach (fr, fatreads) {
		int total = 0;
		// count -> set
		vector<pair<int,int> >  counts;
		foreach (c, fr->clusters) 
			if (clusters[c->first].covered) {
				// count number of reads this set covers
				int cnt = 0;
				foreach (pos, clusters[c->first].fatreads) 
					foreach (fx, pos->second) // in fatreads
						cnt += fatreads[*fx].reads.size();
				total += cnt;
				counts.push_back(make_pair(cnt, c->first));
			}

		sort(counts.begin(), counts.end());
		int remain = total;
		for (int i = 0; i < counts.size(); i++) {
			int r = (counts[i].first * fr->reads.size()) / total;
			fr->solution[ counts[i].second ] += r;
			remain -= r;
		}
		fr->solution[ counts[counts.size() - 1].second ] += remain;
	}
}

/*****************************************************************************************/
/**                                 cplex solver                                        **/
/****************************************************************************************
int smooth (const vector<int> &component) ; 
void cplex_api (IloEnv env, const vector<int> &component, char *model_id) {
	E("\tCPLEX\n");
	// refine the read sets ...
	// [read] -> [[set] -> [var]]
	static vector<map<int, IloNumVar> > read_var;
	if (read_var.size() == 0) {
		read_var.resize(::fatreads.size());
		for (int i = 0; i < ::fatreads.size(); i++) {
			foreach (s, ::fatreads[i].transcripts) 
				read_var[i][s->first] = IloIntVar(env);
		}
	}
	
	//E("Set size %d\n", component.size());
	//fflush(stdout);
	IloModel model(env, model_id);

	// Objective: MINIMIZE SUM(D_transcript)
	IloNumVarArray d(env, component.size(), 0, IloInfinity);
	model.add(IloMinimize(env, IloSum(d)));

	/*int rc = 0;
	foreach (c, component) {
		foreach (p, transcripts[*c].fatreads)
			rc += p->second.size();
	}
	E("Here %d fatreads\n", rc);****-

	// Bounds: D_transcript
	foreachidx (c, ci, component) {
		//IloExprArray nr(env);
		vector<IloExpr> nr;
	
		// get NR
		IloExpr nrsum(env);
		foreach (pos, ::transcripts[*c].fatreads) { 
			nr.push_back(IloExpr(env));
			foreach (r, pos->second)
				nr.back() += read_var[*r][*c]; 
			nrsum += nr.back();
		}

		nrsum /= double(::transcripts[*c].length);
		foreachidx (pos, posi, ::transcripts[*c].fatreads) {
			model.add(nrsum + d[ci] - nr[posi] >= 0);
			model.add(nrsum - d[ci] - nr[posi] <= 0);
		}
		nrsum.end();
		foreach (n, nr)
			n->end();
	}
	
	// Constraint: SUM(Rij) = |R|
	set<int> component_fatreads;
	foreach (c, component) 
		foreach (pos, ::transcripts[*c].fatreads) 
			foreach (r, pos->second)
				component_fatreads.insert(*r);
	foreach (r, component_fatreads) { 
		IloExpr expr(env);
		foreach (s, read_var[*r])
			expr += s->second;
		model.add(expr == ::fatreads[*r].reads.size());
		expr.end();
	}
	component_fatreads.clear();

	//E("Starting CPLEX ...\n");
	IloCplex cplex(model);
	cplex.setParam(IloCplex::Threads, 3);
	cplex.setParam(IloCplex::TiLim, 120);
	cplex.exportModel(model_id);
	cplex.setOut(env.getNullStream());
	bool ok = 0;
	try {
		cplex.solve();
		E("\tSet size %d, cplex status: %d, opt. value %.3lf\n", component.size(), cplex.getCplexStatus(), cplex.getObjValue());
		ok = 1;
	}
	catch (IloException &ex) {
		E("\tException! %s\n", ex.getMessage());
		smooth(component);
	}

	if (ok) foreach (r, component_fatreads) { 
		int sol = 0;
		foreach (s, read_var[*r]) {
			int res = cplex.getValue(s->second);
			if (res)
				::fatreads[*r].solution[s->first] = res;
			sol += res;
		}
		assert(sol == ::fatreads[*r].reads.size());
	}

	E("\tCPLEX %'.2lfM --> ", env.getTotalMemoryUsage()/(1024*1024.0));

	cplex.clearModel();
	cplex.end();
	d.endElements();
	d.end();
	IloExtractableArray del(env);
	for (IloModel::Iterator i(model); i.ok(); ++i)
		del.add(*i);
	del.add(model);
	del.endElements();
	del.end();
	
	E("CPLEX %'.2lfM", env.getTotalMemoryUsage()/(1024*1024.0));
	E("\n");

}

int paramX(char *pre, ...) {
	

}
/*
void cplex_api_file (IloEnv env, const vector<int> &component, char *model_id) {
	E("\tCPLEX\n");
	// refine the read sets ...
	// [read] -> [[set] -> [var]]
	static vector<map<int, IloNumVar> > read_var;
	if (read_var.size() == 0) {
		read_var.resize(::fatreads.size());
		for (int i = 0; i < ::fatreads.size(); i++) {
			foreach (s, ::fatreads[i].transcripts) 
				read_var[i][s->first] = IloIntVar(env);
		}
	}
	
	

	// Objective: MINIMIZE SUM(D_transcript)
	fprintf(lp, "MINIMIZE\n");

	string s = "";

	("D_%d", component[0]);
	for (int i = 0c, component)
		LX(" + D_%d\n", *c);
	fprintf("\n");
	IloNumVarArray d(env, component.size(), 0, IloInfinity);
	model.add(IloMinimize(env, IloSum(d)));

	int rc = 0;
	foreach (c, component) {
		foreach (p, transcripts[*c].fatreads)
			rc += p->second.size();
	}
	E("Here %d fatreads\n", rc);

	// Bounds: D_transcript
	foreachidx (c, ci, component) {
		//IloExprArray nr(env);
		vector<IloExpr> nr;
	
		// get NR
		IloExpr nrsum(env);
		foreach (pos, ::transcripts[*c].fatreads) { 
			nr.push_back(IloExpr(env));
			foreach (r, pos->second)
				nr.back() += read_var[*r][*c]; 
			nrsum += nr.back();
		}

		nrsum /= double(::transcripts[*c].length);
		foreachidx (pos, posi, ::transcripts[*c].fatreads) {
			model.add(nrsum + d[ci] - nr[posi] >= 0);
			model.add(nrsum - d[ci] - nr[posi] <= 0);
		}
		nrsum.end();
		foreach (n, nr)
			n->end();
	}
	
	// Constraint: SUM(Rij) = |R|
	set<int> component_fatreads;
	foreach (c, component) 
		foreach (pos, ::transcripts[*c].fatreads) 
			foreach (r, pos->second)
				component_fatreads.insert(*r);
	foreach (r, component_fatreads) { 
		IloExpr expr(env);
		foreach (s, read_var[*r])
			expr += s->second;
		model.add(expr == ::fatreads[*r].reads.size());
		expr.end();
	}
	component_fatreads.clear();

	//E("Starting CPLEX ...\n");
	IloCplex cplex(model);
	cplex.setParam(IloCplex::Threads, 3);
	cplex.setParam(IloCplex::TiLim, 120);
	cplex.exportModel(model_id);
	cplex.setOut(env.getNullStream());
	bool ok = 0;
	try {
		cplex.solve();
		E("\tSet size %d, cplex status: %d, opt. value %.3lf\n", component.size(), cplex.getCplexStatus(), cplex.getObjValue());
		ok = 1;
	}
	catch (IloException &ex) {
		E("\tException! %s\n", ex.getMessage());
		smooth(component);
	}

	if (ok) foreach (r, component_fatreads) { 
		int sol = 0;
		foreach (s, read_var[*r]) {
			int res = cplex.getValue(s->second);
			if (res)
				::fatreads[*r].solution[s->first] = res;
			sol += res;
		}
		assert(sol == ::fatreads[*r].reads.size());
	}

	E("\tCPLEX %'.2lfM --> ", env.getTotalMemoryUsage()/(1024*1024.0));

	cplex.clearModel();
	cplex.end();
	d.endElements();
	d.end();
	IloExtractableArray del(env);
	for (IloModel::Iterator i(model); i.ok(); ++i)
		del.add(*i);
	del.add(model);
	del.endElements();
	del.end();
	
	E("CPLEX %'.2lfM", env.getTotalMemoryUsage()/(1024*1024.0));
	E("\n");

}

*/

// DANGEROUS ZONE
// I TAKE NO RESPONSIBILITY
// FOR THE CODE BELOW
// (in fact, i am not quite sure what it exactly does... kinda magic ...)

typedef struct { 
	double avg; // average value
	int maxpos; // maxpeak
	map<int, int> poscnt; // graph, X->Y for non-zeros (most of the time)
} CLUSTER_STAT;
#define MIN3(a,b,c)  (min(a,min(b,c)))

/* component -- vector id cluster ids to be smoothed */
int smooth (const vector<int> &component) {
	// TODO: This dares to be UGLY!
	static vector<CLUSTER_STAT> levels;
	if (!levels.size()) 
		levels.resize(clusters.size());

	L("\tHeuristics!\n");
	// total sum of peaks; red is reduced by heuristics
	int sum = 0, red = 0;
	foreach (c, component) {
		// find average and max peak for this cluster
		int l = -1; // max peak
		levels[*c].avg = 0;
		foreach (p, clusters[*c].fatreads) { // p: pos -> [fatread...]
			// sum of weights
			int s = 0; 
			foreach (r, p->second) // r: fatread
				s += fatreads[*r].solution[*c];
			// new max peak?
			if (s > l) {
				l = s;
				levels[*c].maxpos = p->first;
			}
			// value on this position
			// we only consider non-zero positions here
			levels[*c].poscnt[p->first] = s;
			levels[*c].avg += s;
		//	E("(%d -> %d) ", p->first, s);
		}
		levels[*c].avg /= clusters[*c].partial->length;
		sum += l;
	}

	// do it! limit the number of iterations as well
	bool alive = 0;
	int iters = 0;
	do {
		iters++;
		alive = 0;
		// iteration limit
		if (iters > 5000) 
			break;

		foreach (c, component) {
			// <F>rom  component
			int F = *c;
			// its peak
			int c_pos = levels[F].maxpos;
			// weird case? if so, skip! YES this happens a lot
			if (levels[F].poscnt[c_pos] == 0) 
				continue;
			if (levels[F].poscnt[c_pos] <= int(levels[F].avg)) // already avearage?
				continue;

			// now relocate fatreads at the peak position
			foreach (fr, clusters[F].fatreads[c_pos]) if (fatreads[*fr].solution[F]) { // find max, opt it
				// where can we move it?
				foreach (nt, fatreads[*fr].clusters) /* cluster -> pos */ if (clusters[nt->first].covered && nt->first != F) {
					// move <T>o T
					int T = nt->first;
					int t_pos = nt->second;

					// ok; can we move it? any space left?
					if ( fatreads[*fr].solution[F] > 0 && levels[T].poscnt[t_pos] < levels[T].poscnt[levels[T].maxpos] ) {
						// move something to there

						// how many reads should we relocate?
						int how_much = MIN3(
								levels[T].poscnt[levels[T].maxpos] - levels[T].poscnt[t_pos], // destination cap
								fatreads[*fr].solution[F], // source cap
							 	levels[F].poscnt[c_pos] - int(levels[F].avg) // do not go below the source average though  
						);

						//	update results, T part
						fatreads[*fr].solution[T] += how_much;
						levels[T].poscnt[t_pos] += how_much;
						levels[T].avg += double(how_much) / clusters[T].partial->length;
						if (levels[T].poscnt[t_pos] > levels[T].poscnt[levels[T].maxpos]) // update max peak of T
							levels[T].maxpos = c_pos;

						//	update results, F part
						fatreads[*fr].solution[F] -= how_much;
						levels[F].poscnt[c_pos] -= how_much;
						levels[F].avg -= double(how_much) / clusters[F].partial->length;

						if (c_pos == levels[F].maxpos) { // should we update max level of F?
							// ol is old maximum
							int ol = levels[F].poscnt[c_pos] + how_much, l = -1;
							foreach (p, levels[F].poscnt) 
								if (p->second > l) {
									l = p->second;
									levels[F].maxpos = p->first;
								}
							// updated! adjust reduced value
							// TODO maybe this is not correct. who cares anyway?
							red += ol - l;
						}

						// done with F
						if (levels[F].poscnt[c_pos] == 0) 
							goto end;
						// again done with F (this can happen because of rounding)
						if (levels[F].poscnt[c_pos] <= int(levels[F].avg))
							goto end;

						// let's go again!
						alive = 1;
					}
				}
			}
			end:; // goto ftw!
		}
	} while (alive);

//	L("Set sz %'d, iters %'d, ini %'d, red %'d\n", component.size(), iters, sum, red);
}

// DFS from cluster i; n_sets is component big-bang
void connected_component_dfs (int i, vector<int> &n_sets, vector<char> &visited) {
	// add this cluster to the component
	n_sets.push_back(i);
	visited[i] = 1;

	foreach (pos, clusters[i].fatreads) foreach (r, pos->second) { // foreach fatread
		foreach (s, fatreads[*r].clusters) 
			if (!visited[s->first] && clusters[s->first].covered) 
				connected_component_dfs(s->first, n_sets, visited);
	}
}

void connected_components (void) {
//	IloEnv env;

	vector<char> visited(clusters.size(), 0);
	int cc_count = 0, // how many components?
		 cc_sets = 0;  // how many sets covered by components?

	foreach (r, fatreads)
		foreach (s, r->clusters)
			if (!visited[s->first] && clusters[s->first].covered) {	
				vector<int> component;
				connected_component_dfs(s->first, component, visited);

				// calculate number of fatreads in this component
				int rc = 0; 
				foreach (c, component) {
					foreach (p, clusters[*c].fatreads)
						rc += p->second.size();
				}
			
				if (component.size() > 1)
					L("\tBig component: %d transcripts, %d fatreads\n", component.size(), rc);
//				char c[200];
//				sprintf(c, "model_%d.lp", cc_count);
//				cplex_api(env, component, c);
//

				// done, whoa whoa whoa!
				if (component.size() < 2) 
					;
			//	else if (component.size() < 10)
			//		cplex_api(env, component, c);
				else
					smooth(component);
	
				// some stat update
				cc_count++;
				cc_sets += component.size();
			}		

	E("\tCover count %'d [%'d]\n", cc_count, cc_sets);
//	env.end();
}

void update_solution (map<string, struct read> &reads) {
	int discarded = 0;
	foreach (fi, fatreads) {
		int i = 0;
		foreach (si, fi->solution) {
			// debug assertion
			assert(i+si->second <= fi->reads.size() );

			// update reads
			for (int c = 0; c < si->second; c++) {
				struct read *r = fi->reads[i + c];
				cluster &c = clusters[si->first];
				int p_start = fi->clusters[si->first];

				// find the read
				// TODO what if multiple mappings on same pos?
				read::read_key lhs;
				read::read_entry rhs;
				foreach (ei, r->entries)
					if (ei->second.partial == c.partial && ei->second.partial_start == p_start) {
						lhs = ei->first;
						rhs = ei->second;
						break;
					}
				r->entries.clear();
				r->entries[lhs] = rhs;
			}
			i += si->second;
		}
		int no = fi->reads.size() - i;
		assert(no == 0 || no == fi->reads.size());
		for (; i < fi->reads.size(); i++) // remove those reads
			fi->reads[i]->entries.clear();
		discarded += no;
	}
	E("\t%'d discarded reads\n", discarded);
//	fclose(fo);
}

void do_orman (const genome_annotation &ga, const set<partial_transcript> &partials, map<string, struct read> &reads, int read_length) {
	E("Initializing ORMAN; read length is %d ...\n", read_length);
	initialize_structures(partials, reads);
	E("done in %d seconds!\n", zaman_last());

	E("Set cover ...\n");
	int n = set_cover(read_length);
	E("\t%'d out of %'d clusters selected\ndone in %d seconds!\n", n, clusters.size(), zaman_last());

	E("Probabilistic assign ...\n");
	probabilistic_assign();
	E("done in %d seconds!\n", zaman_last());

	E("Connected component decomposition ...\n");
	connected_components();
	E("done in %d seconds!\n", zaman_last());

	E("Updating solution ... \n");
	update_solution(reads);
	E("done in %d seconds!\n", zaman_last());
}

