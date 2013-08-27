/// 786

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
std::string numtostr(int n);

extern bool retain_junk;

#endif // COMMON_H__

