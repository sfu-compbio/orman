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

///
/// node_t<T> definitions
//

template<typename T>
node_t<T>::node_t () {
}

template<typename T>
node_t<T>::node_t (const interval &_k, const T &_v, char c) : 
	key(_k), value(_v), l(0), r(0), color(c) 
{
	max = _k.second;
}

template<typename T>
node_t<T>::~node_t () {
//	delete key; key = 0;
//	delete val; val = 0;

	delete l;
	delete r;
}

template<typename T>
void node_t<T>::setmax () {
	max = key.second;
	if (l && l->max > max) max = l->max;
	if (r && r->max > max) max = r->max;
}

///
/// interval_tree<T> definitions
///

template<typename T>
void interval_tree<T>::color_flip (node_t<T> *n) {
	n->color = !n->color;
	n->l->color = !n->l->color;
	n->r->color = !n->r->color;
}

template<typename T>
node_t<T> *interval_tree<T>::rotate_left (node_t<T> *n) {
	node_t<T> *x = n->r;
	n->r = x->l;
	x->l = n;
	x->color = n->color;
	n->color = RED;

	x->setmax();
	n->setmax();

	return x;
}

template<typename T>
node_t<T> *interval_tree<T>::rotate_right (node_t<T> *n) {
	node_t<T> *x = n->l;
	n->l = x->r;
	x->r = n;
	x->color = n->color;
	n->color = RED;

	x->setmax();
	n->setmax();

	return x;
}

template<typename T>
node_t<T> *interval_tree<T>::insert (node_t<T> *n, const interval &k, const T &v) {
	if (!n)
		return new node_t<T>(k, v, RED);

	if ((n->l && n->l->color == RED) && (n->r && n->r->color == RED))
		color_flip(n);

	if (k.first < n->key.first) 
		n->l = insert(n->l, k, v);
	else if (k.first == n->key.first && k.second == n->key.second && v == n->value)
		;
	else 
		n->r = insert(n->r, k, v);

	if (n->r && n->r->color == RED)
		n = rotate_left(n);
	if ((n->l && n->l->color == RED) && (n->l && n->l->l && n->l->l->color == RED))
		n = rotate_right(n);

	n->setmax();
	return n;
}

template<typename T>
void interval_tree<T>::enumerate (interval_range p, std::vector<T> &v, node_t<T> *n) {
	if (!n || p > n->max)
		return;
	enumerate(p, v, n->l);
	if (p >= n->key.first && p <= n->key.second)
		v.push_back(n->value);
	if (p < n->key.first)
		return;
	enumerate(p, v, n->r);
}

template<typename T>
void interval_tree<T>::enumerate (const interval &p, std::vector<T> &v, node_t<T> *n) {
	if (!n || p.first > n->max)
		return;
	enumerate(p, v, n->l);
	if (n->key.first <= p.second && n->key.second >= p.first)
		v.push_back(n->value);
	if (p.second < n->key.first)
		return;
	enumerate(p, v, n->r);
}

template<typename T>
void interval_tree<T>::insert (const interval &p, const T &v) {
	root = insert(root, p, v);
	root->color = BLACK;
}

template<typename T>
void interval_tree<T>::enumerate (interval_range p, std::vector<T> &v) {
	enumerate(p, v, root);
}

template<typename T>
void interval_tree<T>::enumerate (const interval &p, std::vector<T> &v) {
	enumerate(p, v, root);
}

template<typename T>
interval_tree<T>::interval_tree () : 
	root(0) 
{
}

template<typename T>
interval_tree<T>::~interval_tree () { 
	delete root; 
}

#endif
