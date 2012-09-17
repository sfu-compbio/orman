/// 786

#include "common.h"
#include <boost/heap/fibonacci_heap.hpp>
#include <ilcplex/ilocplex.h>

using namespace std;
using namespace boost::heap;


FILE *flog;
#define LE(c,...) fprintf(flog,c,##__VA_ARGS__)

struct transcript_t {
	// ID = Position in ::fatreads array
	int id;
	string gene, signature;
	// Length of transcript, len=0 is intronic
	int length, weight;
	// Links to the fatreads, [position] -> [<fatread_index, ...>]
	// TODO vector ili set?
	map<int, set<int> > fatreads;

public: 
	transcript_t (void) {}
	transcript_t (int i, const char *G, const char *S, int l, int w):
		id(i), gene(G), signature(S), length(l), weight(w) {}
};
vector<transcript_t> transcripts;

struct fatread_t {
	// IDs of the reads
	vector<string> reads;
	// Transcript ID -> position
	// Transcript ID SHOULD be unique
	map<int, int> transcripts;

	// TrID -> HowMuch
	map<int, int> solution;

public:
	fatread_t (void): reads(0) {} 
};
vector<fatread_t> fatreads;

/*****************************************************************************************/
/**                                    reading                                          **/
/*****************************************************************************************/
void read_input (const char *path) {
	char buffep[MAX_BUFFER],
		  buffer[MAX_BUFFER];
	FILE *fi = fopen(path, "r");

	// read fatreads
	map<string, map<int, int> > reads_index;
	while (fscanf(fi, "%s", buffer)) {
		if (buffer[0] == '!')
			break;

		map<int, int> tx;

		int ns; fscanf(fi, "%d", &ns);
		for (int i = 0; i < ns; i++) {
			int set, pos; fscanf(fi, "%d %d", &set, &pos);
			if (set != -1) 
				tx[set] = pos;
		}

		auto it = reads_index.insert(make_pair(string(buffer), tx)).first;
	}
	E("\tRead %'d distinct reads\n", reads_index.size());

	int ns, si, sl, sle, swe;
	fscanf(fi, "%d", &ns);
	::transcripts.resize(ns); 
	while (fscanf(fi, "%d %d %s %s %d %d", &si, &sl, buffep, buffer, &sle, &swe) != EOF)
		::transcripts[si] = transcript_t(si, buffep, buffer, sle, swe);
	E("\tRead %'d transcripts\n", ::transcripts.size());

	map<string, int> fatread_index;
	foreach (r, reads_index) {
		if (r->second.size() < 1 || r->second.size() > 100)
			continue;

		string fatread_signature = "";
		foreach (pp, r->second) {
			sprintf(buffep, "%d+%d;", pp->first, pp->second);
			fatread_signature += string(buffep);
		}

		auto ft = fatread_index.find(fatread_signature);
		if (ft == fatread_index.end()) {
			ft = fatread_index.insert(make_pair(fatread_signature, fatread_index.size())).first;
			if (::fatreads.size() <= ft->second) 
				::fatreads.resize(ft->second + 1);
			foreach (pp, r->second) {
				::fatreads[ft->second].transcripts[pp->first] = pp->second;
				::fatreads[ft->second].solution[pp->first] = 0;
			}
		}
		::fatreads[ft->second].reads.push_back(r->first);
	}
	E("\tFormed %'d fatreads\n", ::fatreads.size());

	foreachidx (r, ri, ::fatreads) {
		foreach (t, r->transcripts) 
			::transcripts[t->first].fatreads[t->second].insert(ri);
	}

	fclose(fi);
}



/*****************************************************************************************/
/**                                    set cover                                        **/
/*****************************************************************************************/
struct heap_t {
	transcript_t *t;

public:
	int cardinality;
	double weight;

public:
	heap_t (void): 
		t(0) {}
	heap_t (transcript_t *t, int c, double w): 
		t(t), cardinality(c), weight(w) {}

public:
	double value (void) const {
		return cardinality / weight;
	}
	bool operator< (const heap_t &b) const {
		return value() < b.value();
	}
};

vector<char> in_cover;
int cover (void) {
	fibonacci_heap<heap_t> heap;
	vector<fibonacci_heap<heap_t>::handle_type> heap_handles(::transcripts.size());
	for (int i = 0; i < ::transcripts.size(); i++) {
		int covers = 0;
		// count number of reads this set covers
		foreach (pos, ::transcripts[i].fatreads) 
			foreach (r, pos->second)
				covers += ::fatreads[*r].reads.size();
		heap_handles[i] = heap.push(heap_t(&::transcripts[i], covers, ::transcripts[i].weight));
	}

	vector<char> read_visited(::fatreads.size(), 0);
	in_cover = vector<char>(::transcripts.size(), 0);
	int mincover = 0;
	while (!heap.empty()) {
		heap_t h = heap.top(); heap.pop();

		if (h.cardinality == 0)
			break;
		mincover++;

		foreach (pos, h.t->fatreads) foreach (r, pos->second) 
			if (!read_visited[*r]) {
				::fatreads[*r].solution[ h.t->id ] += ::fatreads[*r].reads.size(); // assign ALL to this set

				foreach (s, ::fatreads[*r].transcripts) if (s->first != h.t->id) {
					(*heap_handles[s->first]).cardinality -= ::fatreads[*r].reads.size();
					heap.decrease(heap_handles[s->first]);
				}
				read_visited[*r] = 1;
			}
		in_cover[h.t->id] = 1;
	//	L("%10d[->%20s,%20s] %10d =[%10.0lf]\n", h.t->id, h.t->gene.c_str(), h.t->signature.c_str(), _tx, h.value());
	}

	int not_covered = 0, ncr = 0;
	foreach (fr, ::fatreads) {
		char cv = 0;
		foreach (t, fr->transcripts)
			cv |= in_cover[t->first];;
		if (!cv) {
			not_covered++;
			ncr += fr->reads.size();
		}
	}
	E("\t%'d non-covered fatreads (%'d reads)\n", not_covered, ncr);

	return mincover;
}

/*****************************************************************************************/
/**                                 cplex solver                                        **/
/*****************************************************************************************/
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
	E("Here %d fatreads\n", rc);*/

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
	cplex.setParam(IloCplex::TiLim, 5);
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

typedef struct { 
	double avg;
	int maxpos; 
	map<int, int> poscnt; 
} STAT;
int smooth (const vector<int> &component) {
	// [[pos]->[exp]] [maxpos]
	
	static vector<STAT> levels;
	if (!levels.size()) 
		levels.resize(transcripts.size());

	E("\tHeuristics!\n");
	LE(">component %d\n", component.size());
	foreach(c, component)
		LE("\t>transcript %s -> %s\n", transcripts[*c].gene.c_str(), transcripts[*c].signature.c_str());

	int sum = 0, red = 0;
	foreach (c, component) {
		int l = -1;
		levels[*c].avg = 0;
		foreach (p, transcripts[*c].fatreads) { // p: pos -> [fatread...]
			int s = 0;
			foreach (r, p->second) { // r: fatread
				assert( fatreads[*r].solution.find(*c) != fatreads[*r].solution.end() );
				s += fatreads[*r].solution[*c];
			}
			if (s > l) {
				l = s;
				levels[*c].maxpos = p->first;
			}
			levels[*c].poscnt[p->first] = s;
			levels[*c].avg += s;
		}
		levels[*c].avg /= transcripts[*c].length;
		sum += l;
	}
//	E("sz=%d d=%'d\n",component.size(), sum);

	bool alive = 0;
	int iters = 0;
	do {
		iters++;
		alive = 0;
		if (iters > 5000) break;

		foreach (c, component) {

			int F = *c;

			int c_pos = levels[F].maxpos;
			if (levels[F].poscnt[c_pos] == 0) continue;

			// take first one
//			assert ( transcripts[*c].fatreads.find(c_pos) != transcripts[*c].fatreads.end() );
//			E(">tr.sz=%d\n",transcripts[*c].fatreads[c_pos].size());

			foreach (fr, transcripts[*c].fatreads[c_pos]) if (fatreads[*fr].solution[*c]) { // find max, opt it
//				E("Trying to relocate %d ~~%d\n", *fr, fatreads[*fr].solution[*c]);
				// where can we move it?

				foreach (nt, fatreads[*fr].transcripts) if (in_cover[nt->first] && nt->first != *c) {
					int T = nt->first;
					int t_pos = nt->second;

					if ( levels[T].poscnt[t_pos] < levels[T].poscnt[levels[T].maxpos] ) {
						// move something to there

						/*
						 *
						 * MOVE up to current average
						 * UPDATE AVERAGE!
						 *
						 */

//						assert( levels[F].poscnt.find(c_pos) != levels[F].poscnt.end() );
//						assert( levels[T].poscnt.find(c_pos) != levels[T].poscnt.end() );
						int how_much = min(levels[T].poscnt[levels[T].maxpos] - levels[T].poscnt[t_pos], // possible
												 fatreads[*fr].solution[*c] - int(levels[F].avg) /* levels[F].poscnt[c_pos] */); // remaining
						
//						if(how_much<=0) {
//							E( "M( %d - %d, %d ) %d[%d] %d | %d %d\n", levels[T].poscnt[levels[T].maxpos], levels[T].poscnt[c_pos], levels[F].poscnt[c_pos], F, *c, T, *c, c_pos );
//							exit(1);
//						}
//						E("Moding %d from %d to %d, pos %d->%d\n", how_much, F, T, c_pos, t_pos);						
//						E("O %d %d\n", fatreads[*fr].solution[F], fatreads[*fr].solution[T]);
						fatreads[*fr].solution[F] -= how_much;
						levels[F].avg -= double(how_much) / transcripts[F].length;
						fatreads[*fr].solution[T] += how_much;
						levels[T].avg += double(how_much) / transcripts[T].length;
//						E("N %d %d\n", fatreads[*fr].solution[F], fatreads[*fr].solution[T]);

//						assert( levels[T].poscnt.find(c_pos) != levels[T].poscnt.end() );
						levels[T].poscnt[t_pos] += how_much;
						if (levels[T].poscnt[t_pos] > levels[T].poscnt[levels[T].maxpos]) // update TO max level
							levels[T].maxpos = c_pos;


//						assert( levels[F].poscnt.find(c_pos) != levels[F].poscnt.end() );
						levels[F].poscnt[c_pos] -= how_much;
//						levels[F].sumcnt -= how_much;

						if (c_pos == levels[F].maxpos) {
							int ol = levels[F].poscnt[c_pos] + how_much, l = -1;
							foreach (p, levels[F].poscnt) {
								if (p->second > l) {
									l = p->second;
									levels[F].maxpos = p->first;
								}
							}
							//						assert(ol >= l);
							red += ol - l;
						}

						if (levels[F].poscnt[c_pos] == 0) goto end;

						alive = 1;
					}
				}
			}
end:;
		}

	} while (alive);

	LE("\t>iters %d, reduced %d\n", iters, red);

//	L("Set sz %'d, iters %'d, ini %'d, red %'d\n", component.size(), iters, sum, red);

	//	L("ROW %'d COLSUM %'d\n", nr, nc);
}

void connected_component_dfs (int i, vector<int> &n_sets, vector<char> &visited) {
	n_sets.push_back(i);
	visited[i] = 1;


	foreach (pos, ::transcripts[i].fatreads) foreach (r, pos->second) {
		foreach (s, ::fatreads[*r].transcripts) 
			if (!visited[s->first] && in_cover[s->first]) {
				
			//	L("Node %d: %s+%s\n", i, transcripts[i].gene.c_str(), transcripts[i].signature.c_str());
			//	L("\tvia "); 

				connected_component_dfs(s->first, n_sets, visited); 
			}
	}
}

void connected_components (void) {
	vector<char> visited(::transcripts.size(), 0);

	IloEnv env;

	int cc_count = 0,
		 cc_sets = 0;

	E("\n");

	foreach (r, ::fatreads)
		foreach (s, r->transcripts)
			if (!visited[s->first] && in_cover[s->first]) {	
				vector<int> component;
				connected_component_dfs(s->first, component, visited);

				int rc = 0;
				foreach (c, component) {
					foreach (p, transcripts[*c].fatreads)
						rc += p->second.size();
				}
			
				if (component.size() > 1)
					E("Component %d transcripts, %d fatreads\n", component.size(), rc);
				char c[200];
				sprintf(c, "model_%d.lp", cc_count);
//				cplex_api(env, component, c);
//

				if (component.size() < 2) 
					;
				else if (component.size() < 10)
					cplex_api(env, component, c);
				else
					smooth(component);
	
				cc_count++;
				cc_sets += component.size();
			}		

	E("\n\tCover count %'d [%'d]\n", cc_count, cc_sets);
	env.end();
}

void print_solution () {
//	FILE *fo = fopen(out.c_str(), "w");
	int no = 0;
	foreach (f, fatreads) {
		int i = 0;
		foreach (s, f->solution) {
			//if (!s->second)
				
			for (int c = 0; c < s->second; c++)
				L("%s %d %d\n", 
						f->reads[i + c].c_str(), 
						s->first, 
						f->transcripts[s->first]
				);
			i += s->second;
		}
		no += f->reads.size() - i;
		for (; i < f->reads.size(); i++)
			L("%s -1\n", f->reads[i].c_str()); 
	}
	E("\t%'d discarded reads\n", no);
//	fclose(fo);
}

int main (int argc, char **argv) {
	setlocale(LC_ALL, "");
	zaman_last();

	flog = fopen("stage2.log", "w");

	char buffer[MAX_BUFFER];
	E("Reading stage1 output from %s ... \n", realpath(argv[1], buffer));
	read_input(argv[1]);
	E("done in %d seconds!\n", zaman_last());

	E("Set cover ...\n");
	int n = cover();
	E("\t%'d transcripts found\ndone in %d seconds!\n", n, zaman_last());

	E("Connected component decomposition ...        ");
	connected_components();
	E("done in %d seconds!\n", zaman_last());

	E("Printing solution ... \n");
	print_solution();
	E("done in %d seconds!\n", zaman_last());

	fclose(flog);

	return 0;
}

