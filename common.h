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

#ifndef COMMON_H__
#define COMMON_H__

#include <ctime>
#include <set>
#include <map>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <iostream>
#include <cstdio>
#include <vector>
#include <locale>
#include <cstdlib>
#include <cctype>
#include <ctime>
#include <string>
#include <cstring>
#include <algorithm>

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

#define KB (1024)
#define MB (KB*KB)
#define GB (KB*MB)

#define W()      getchar()
#define L(c,...) fprintf(stdout,c,##__VA_ARGS__)
#define E(c,...) fprintf(stderr,c,##__VA_ARGS__)

#define VALID(c) ((c)!='N')

#define foreach(i,c) \
  for (typeof((c).begin()) i = (c).begin(); i != (c).end(); ++i)
#define foreachidx(i,I,c) \
	int I = 0; for (typeof((c).begin()) i = (c).begin(); i != (c).end(); ++i, I++) 

const int MAX_BUFFER = 25000;
uint64_t zaman_start();
uint64_t zaman_last();
string numtostr(int n);

extern int optThreads;

#endif // COMMON_H__

