/// 786

#ifndef RESCUE_H__
#define RESCUE_H__

#include "common.h"
#include "partial.h"
#include "annotation.h"

void do_single (const genome_annotation &ga, map<std::string, struct read> &reads);
void do_rescue (const genome_annotation &ga, map<std::string, struct read> &reads);

#endif // COMMON_H__

