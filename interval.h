/// 786

#ifndef INTERVAL_H__
#define INTERVAL_H__

#include <inttypes.h>
#include <vector>
#include <utility>

typedef uint32_t interval_range;
typedef std::pair<interval_range, interval_range> interval;

const char RED   = 0;
const char BLACK = 1;

template<typename T>
struct node_t {
	interval key;
	T value;

	interval_range max;

	char color;
	node_t<T> *l, *r;

public:
	node_t ();
	node_t (const interval &_k, const T &_v, char c);
	~node_t ();

	void setmax ();
};

template<typename T>
class interval_tree {
	node_t<T> *root;
	
private:
	void color_flip (node_t<T> *n);
	node_t<T> *rotate_left (node_t<T> *n);
	node_t<T> *rotate_right (node_t<T> *n);
	node_t<T> *insert (node_t<T> *n, const interval &k, const T &v);
	void enumerate (interval_range p,  std::vector<T> &v, node_t<T> *n);
	void enumerate (const interval &p, std::vector<T> &v, node_t<T> *n);
	
public:	
	interval_tree ();
	~interval_tree ();

	void insert (const interval &k, const T &v);
	void enumerate (interval_range p,  std::vector<T> &result);
	void enumerate (const interval &p, std::vector<T> &result);
};

#include "interval.tc"

#endif
