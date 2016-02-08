/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_READANNO_H_
#define _DD_READANNO_H_

#include "seq.h"
#include "drompa_gv.h"
#include "dd_gv.h"

void read_graph(Graph *graph, RefGenome *g, int chr, char *name, double min, double max, Graph_Type type);
void read_gene(GeneSet *p, RefGenome *g, int chr, GeneFile_Type gftype);
void read_ARS_OriDB(char *argv, GeneSet *p, int chr);
void read_TER_scer(char *argv, GeneSet *p, int chr);
void read_repeat(DDParam *p, RefGenome *g, int chr);
void read_genefile(DDParam *d, BedFile *p, RefGenome *g, int chr);

#endif /* _DD_READANNO_H_ */
