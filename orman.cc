/// 786

#include "common.h"
#include "orman.h"
#include <boost/heap/fibonacci_heap.hpp>
#include <ilcplex/ilocplex.h>
#include <inttypes.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <pthread.h>

using namespace std;
using namespace boost::heap;

struct cluster {
	int id;
	PT  partial;
	// position -> [<fatreads ...>]
	map<int, set<int> > fatreads;
	vector<int> single; // single coverage
	int  coverage;
	bool covered;

	cluster (void) {}
	cluster (int i, const PT &p) :
		id(i), partial(p), covered(0), coverage(0) {}
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

string print_pt(const PT &p)  {
	char c[500];
	PTsp pt1 = p.first, pt2=p.second;
	sprintf(c, "%s_%s", string(pt1->transcript->gene->name + "." + pt1->signature).c_str(), pt2?string(pt1->transcript->gene->name + "." + pt1->signature).c_str():"");
	return string(c);
}

string print_pt(const cluster *p)  {
	return print_pt(p->partial);
}

/*****************************************************************************************/

int partial_length (const PT &p)  {
	return p.first->length;
}

int partial_weight (const PT &p) {
	if (!p.second)
		return 1000000 + p.first->weight;
	else if (p.first->transcript == p.second->transcript)
		return p.first->weight + p.second->weight;
	else 
		return p.first->weight + p.second->weight + 100000;
}

void initialize_structures (vector<struct read> &reads) {
	// assign clusters to partials
	set<PT> partials;
	foreach (ri, reads)
		foreach (re, ri->entries) 
		partials.insert(re->partial);
	clusters.resize(partials.size());
	int id = 0;
	map<PT, int> cluster_index;
	foreach (ci, partials) {
		clusters[id] = cluster(id, *ci);
		clusters[id].single.resize(partial_length(*ci));
		for (int k = 0; k < clusters[id].single.size(); k++) {
			clusters[id].single[k] = get_single_coverage(*ci, k);
			clusters[id].coverage += clusters[id].single[k];
		}
		cluster_index[*ci] = id;
		id++;
	}
	partials.clear();

	// initialize fatreads
	char buffer[MAX_BUFFER];
	map<string, int> fatread_index;
	foreach (ri, reads) {
		struct read &r = *ri;
		if (r.entries.size() < 1)
			continue;

		// prepare and sort signature!
		// SIG: <PT+PT1_loc>
		map<PT, int> sig;
		foreach (ei, r.entries) 
			sig.insert(make_pair(ei->partial, ei->partial_start.first));
		// form signature
		string signature = "";
		foreach (si, sig) {
			sprintf(buffer, "%llu+%d;", cluster_index[si->first], si->second);
			signature += buffer;
		}

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
			clusters[ci->first].fatreads[ci->second].insert(fi);
		}
	}

	// initialize coverage
	for (int i = 0; i < clusters.size(); i++) {
		foreach (pos, clusters[i].fatreads) 
			foreach (fr, pos->second) // in fatreads
			clusters[i].coverage += fatreads[*fr].reads.size();
	}

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
		covers = clusters[i].coverage;
		// calculate number of non-covered positions!
		// and update the weight
		double weight = partial_weight(clusters[i].partial);
		int notcovered = 0;
		int prev = 0;
		foreach (fr, clusters[i].fatreads) {
			if (fr->first > prev)
				notcovered += fr->first - prev;
			prev = fr->first + read_length;
		}
		if (prev < partial_length(clusters[i].partial))
			notcovered += partial_length(clusters[i].partial) - prev;
		if (notcovered > 1)
			weight *= (1 + 10.0 * double(notcovered) / partial_length(clusters[i].partial));
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
			foreach (r, pos->second) // in fatreads 
			if (!read_visited[*r]) {
				foreach (s, fatreads[*r].clusters) if (s->first != h.t->id) {
					(*heap_handles[s->first]).cardinality -= fatreads[*r].reads.size();
					(*heap_handles[s->first]).t->coverage -= fatreads[*r].reads.size();
					heap.decrease(heap_handles[s->first]);
				}
				read_visited[*r] = 1;
			}
		h.t->covered = true;
	}

	return clusters.size();
}

void probabilistic_assign (void) {
	// assign farteads to clusters
	// based on cluster coverage
	// higher coverage => more reads
	foreach (fr, fatreads) {
		int total = 0;
		// count -> set
		vector<pair<int,int> >  counts;
		foreach (c, fr->clusters) 
			if (clusters[c->first].covered) {
				int cnt = 0;
				cnt = clusters[c->first].coverage;
				total += cnt;
				counts.push_back(make_pair(cnt, c->first));
			}

		// lo count->hi count
		sort(counts.begin(), counts.end());
		int remain = fr->reads.size();
		for (int i = 0; i < counts.size(); i++) {
			int r = (counts[i].first * fr->reads.size()) / total;
			fr->solution[ counts[i].second ] += r;
			remain -= r;
		}
		fr->solution[ counts[counts.size() - 1].second ] += remain;
	}
}

static pthread_t *threads;
static pthread_mutex_t mtx_access;
static pthread_mutex_t mtx_io;
static int thread_counter = 0;
static vector<vector<int> > components;
static vector<int> ids; 
class raii_lock { // thanks http://stackoverflow.com/questions/14272969/how-to-check-whether-pthread-mutex-t-is-locked-or-not
    pthread_mutex_t &mutex;
public:   
    raii_lock(pthread_mutex_t &m): mutex(m) { pthread_mutex_lock(&mutex); }
    ~raii_lock(void) { pthread_mutex_unlock(&mutex); }
};


typedef struct { 
	double avg; // average value
	int maxpos; // maxpeak
	int minpos;
	vector<int> poscnt; // graph, X->Y for non-zeros (most of the time)
	vector<int> poscnt_; // graph, X->Y for non-zeros (most of the time)
} CLUSTER_STAT;
#define MIN3(a,b,c)  (min(a,min(b,c)))
vector<CLUSTER_STAT> levels;

/* component -- vector id cluster ids to be smoothed */
int heuristics_smooth (const vector<int> &component, int _id, int iter_limit = 5000) {
	pthread_mutex_lock(&mtx_io);
	L("#Heuristics!\n");
	fflush(stdout);
	pthread_mutex_unlock(&mtx_io);
	// total sum of peaks; red is reduced by heuristics
	int sum = 0, red = 0;
	foreach (c, component) {
		// find average and max peak for this cluster
		int l = -1; // max peak
		int L = INT_MAX; // min peak
		levels[*c].avg = 0;
		levels[*c].poscnt.resize(clusters[*c].single.size());
		levels[*c].poscnt_.resize(clusters[*c].single.size());
		for (int i=0;i<clusters[*c].single.size();i++)
	 		levels[*c].poscnt[i] = levels[*c].poscnt_[i] = clusters[*c].single[i];
		foreach (p, clusters[*c].fatreads) { // p: pos -> [fatread...]
			// sum of weights
			int s = 0; 
			foreach (r, p->second) // r: fatread
				s += fatreads[*r].solution[*c];
			// new max peak?
			s += clusters[*c].single[p->first];
			if (s > l) {
				l = s;
				levels[*c].maxpos = p->first;
			}
			if (s < L) {
				L = s;
				levels[*c].minpos = p->first;
			}
			s -= clusters[*c].single[p->first];
			// value on this position
			// we only consider non-zero positions here
			levels[*c].poscnt[p->first] += s;
			levels[*c].poscnt_[p->first] += s;
			//	E("(%d -> %d) ", p->first, s);
		}
		for (int i=0;i<clusters[*c].single.size();i++)
			levels[*c].avg += levels[*c].poscnt[i];
		levels[*c].avg /= partial_length(clusters[*c].partial);
		sum += l;
	}

	// do it! limit the number of iterations as well
	bool alive = 0;
	int iters = 0;
	do {
		iters++;
		alive = 0;
		// iteration limit
		if (iters > iter_limit) 
			break;

		foreach (c, component) {
			// <F>rom  component
			int F = *c;
			// its peak
			int c_pos = levels[F].maxpos;
			// weird case? if so, skip! YES this happens a lot
			if (levels[F].poscnt[c_pos] <= clusters[F].single[c_pos]) 
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
						//ADDED NEW!!!
						how_much = min(how_much, levels[F].poscnt[c_pos] - clusters[F].single[c_pos]);

						//	update results, T part
						fatreads[*fr].solution[T] += how_much;
						levels[T].poscnt[t_pos] += how_much;
						levels[T].avg += double(how_much) / partial_length(clusters[T].partial);
						if (levels[T].poscnt[t_pos] > levels[T].poscnt[levels[T].maxpos]) // update max peak of T
							levels[T].maxpos = c_pos;

						//	update results, F part
						fatreads[*fr].solution[F] -= how_much;
						levels[F].poscnt[c_pos] -= how_much;
						levels[F].avg -= double(how_much) / partial_length(clusters[F].partial);

						if (c_pos == levels[F].maxpos) { // should we update max level of F?
							// ol is old maximum
							int ol = levels[F].poscnt[c_pos] + how_much, l = -1;
							for (int _p = 0; _p < levels[F].poscnt.size(); _p++) 
								if (levels[F].poscnt[_p] > l && levels[F].poscnt[_p] > clusters[F].single[_p]) {
									l = levels[F].poscnt[_p];
									levels[F].maxpos = _p;
								}
							// updated! adjust reduced value
							// TODO maybe this is not correct. who cares anyway?
							red += ol - l;
						}

						// done with F
						if (levels[F].poscnt[c_pos] <= clusters[F].single[c_pos]) 
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
	/*foreach (c, component) {
		L("# %s\n", print_pt(clusters[*c].partial).c_str());
		for (int i=0;i<levels[*c].poscnt_.size();i++) L("(%d,%d) ",i,levels[*c].poscnt_[i]); L("\n");
		for (int i=0;i<levels[*c].poscnt_.size();i++) L("(%d,%d) ",i,levels[*c].poscnt [i]); L("\n");
	}*/
}

void cplex_smooth (const vector<int> &component, int id) {
	IloEnv env;

	// component -> (fatread)
	map<pair<int, int>, IloIntVar> variables;
	for (int i = 0; i < component.size(); i++)
		foreach (fr_pos, clusters[component[i]].fatreads)
			foreach (fr, fr_pos->second) {
				if (variables.find(make_pair(*fr, component[i])) == variables.end())
					variables[make_pair(*fr, component[i])] = IloIntVar(env);
			}

	char model_id[500];
	sprintf(model_id, "orman%06d.lp", id);

	IloModel model(env, model_id);

	// Constraint
	IloNumVarArray d(env, component.size(), 0, IloInfinity);
	model.add(IloMinimize(env, IloSum(d)));
	// -D <= AV - NR <= D
	for (int ci = 0; ci < component.size(); ci++) {
		int c = component[ci];

		vector<IloExpr> nr;
		IloExpr avg(env);
		foreach (pos, clusters[c].fatreads) { 
			nr.push_back(IloExpr(env));
			nr.back() += (double)clusters[c].single[pos->first];
			foreach (r, pos->second)
				nr.back() += variables[make_pair(*r, c)];
			avg += nr.back();
		}

		for (int i = 0; i < clusters[c].single.size(); i++)
			if (clusters[c].fatreads.find(i) == clusters[c].fatreads.end()) 
				avg += (double)clusters[c].single[i];
		avg /= (double)partial_length(clusters[c].partial);
		foreachidx (pos, posi, clusters[c].fatreads) {
			model.add(avg + d[ci] - nr[posi] >= 0);
			model.add(avg - d[ci] - nr[posi] <= 0);
		}
		avg.end();
		foreach (n, nr) n->end();
	}	
	// SUM(Rij) = sizeof(R)
	int fr = 0;
	auto v = variables.begin();
	while (v != variables.end()) {
		IloExpr expr(env);
		int read = v->first.first;
		while (v != variables.end() && v->first.first == read) {
			expr += v->second; 
			v++;
		}
		fr++;
		model.add(expr == fatreads[read].reads.size());
		expr.end();
	}

	{
		raii_lock _l(mtx_io);
		L("Set size %'d, fatreads %'d\n", component.size(), fr);
		fflush(stdout);
	}

	try {
		IloCplex cplex(model);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::TiLim, 120);
//		cplex.setParam(IloCplex::TreLim, 8162);
//		cplex.setParam(IloCplex::WorkMem, 8162);
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
//		cplex.exportModel((string("models/") + model_id).c_str());
		
		cplex.solve();
		
		{
			raii_lock _l(mtx_io);
			L("#CPLEX Status: %d; Optimum value: %.3lf; Memory: %'.2lfM; Instance: %s;\n", cplex.getCplexStatus(), cplex.getObjValue(), env.getTotalMemoryUsage()/(1024*1024.0), model_id);
			fflush(stdout);
		}

		// update solution
		foreach (var, variables) {
			int res = IloRound(cplex.getValue(var->second));
			fatreads[var->first.first].solution[var->first.second] = res;
		}
		// clean
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
	}
	catch (IloException &ex) {
		{
			raii_lock _l(mtx_io);
			string s = ex.getMessage();
			s[s.size() - 1] = 0; // trim \n
			L("#CPLEX Exception: %s; Instance: %s; ", s.c_str(), model_id);
			fflush(stdout);
		}
		ex.end();
		heuristics_smooth(component, id);
	}
	env.end();
}

void *thread (void *arg) {
	while (1) {
		bool survive = true;
		pthread_mutex_lock(&mtx_access);
			if (components.size () == 0) {
				pthread_mutex_unlock(&mtx_access);
				return 0;	
			}
			vector<int> component = components.back();
			int id = ids.back();
			components.pop_back();
			ids.pop_back();
		pthread_mutex_unlock(&mtx_access);

		int rc = 0; 
		foreach (c, component) {
			foreach (p, clusters[*c].fatreads)
				rc += p->second.size();
		}

		if (component.size () <= 1 || rc <= 2) 
			continue;

		cplex_smooth(component, id);
	//	heuristics_smooth(component, id, 1000);
	}
	return 0;
}

// DFS from cluster i; n_sets is component big-bang
void connected_component_dfs (int i, vector<int> &n_sets, vector<char> &visited) {
	n_sets.push_back(i);
	visited[i] = 1;

	foreach (pos, clusters[i].fatreads) foreach (r, pos->second) { // foreach fatread
		foreach (s, fatreads[*r].clusters) 
			if (!visited[s->first] && clusters[s->first].covered) 
				connected_component_dfs(s->first, n_sets, visited);
	}
}

void connected_components (void) {
	vector<char> visited(clusters.size(), 0);
	int cc_count = 0, // how many components?
		 cc_sets = 0; // how many sets covered by components?

	foreach (r, fatreads) foreach (s, r->clusters)
		if (!visited[s->first] && clusters[s->first].covered) {	
			vector<int> component;
			connected_component_dfs(s->first, component, visited);

			components.push_back(component);
			ids.push_back(cc_count++);
			cc_sets += component.size();
		}		

	E("\t%'d components (totaling %'d clusters)\n", cc_count, cc_sets);

	levels.resize(clusters.size());

	// start threading!!
	int thread_count = min((int)sysconf(_SC_NPROCESSORS_ONLN) - 1, 20);
	threads = new pthread_t[thread_count];
	pthread_mutex_init(&mtx_access, 0);
	pthread_mutex_init(&mtx_io, 0);

	for (int i = 0; i < thread_count; i++)
		pthread_create(&threads[i], 0, thread, 0);
	while (1) {
		bool survive = true;
		int size;
		pthread_mutex_lock(&mtx_access);
			if ((size = components.size()) == 0)
				survive = false;
		pthread_mutex_unlock(&mtx_access);
		if (!survive) break;
		
		E("\r\t(%d threads) %'d of %'d completed (%.2lf%%) ...", thread_count, cc_count-size, cc_count, 100.0*(cc_count-size)/double(cc_count));
		sleep(1);
	}
	E("\r\t%'d of %'d completed (%.2lf%%)", cc_count-components.size(), cc_count, 100.0*(cc_count-components.size())/double(cc_count));
	E("\n");
	for (int i = 0; i < thread_count; i++)
		pthread_join(threads[i], 0);

	pthread_mutex_destroy(&mtx_access);
	pthread_mutex_destroy(&mtx_io);
	delete[] threads;
}

void update_solution (vector<struct read> &reads) {
	int discarded = 0;
	foreach (fi, fatreads) {
		int i = 0;
		foreach (si, fi->solution) {
			// debug assertion
			assert(i+si->second <= fi->reads.size() );

			// update reads
			for (int c = 0; c < si->second; c++) {
				struct read *r = fi->reads[i + c];
				cluster &cx = clusters[si->first];
				int p_start = fi->clusters[si->first];

				// find the read
				// TODO what if multiple mappings on same pos?
				read::read_entry rhs;
				foreach (ei, r->entries)
					if (ei->partial == cx.partial && ei->partial_start.first == p_start) {
						rhs = *ei;
						break;
					}
				r->entries.clear();
				r->entries.push_back(rhs);
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

void do_orman (const genome_annotation &ga, vector<struct read> &reads, int read_length) {
	E("Initializing ORMAN; read length is %d ...\n", read_length);
	initialize_structures(reads);
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

