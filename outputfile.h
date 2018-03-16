/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _OUTPUTFILE_H_
#define _OUTPUTFILE_H_

#include "seq.h"
#include "common.h"
#include "drompa_gv.h"

void output_bindata(char *dir, char *prefix, RefGenome *g, TYPE_WIGARRAY *array, char *gtfile, int binsize, int binnum, int chr, PWfile_Type wtype, int showzero);
void make_binary(TYPE_WIGARRAY *array, char *outputfile, int binnum);
void make_bedGraph(RefGenome *g, TYPE_WIGARRAY *array, char *outputfile, char *prefix, int binsize, int binnum, int chr, int showzero);
void convert_bedGraph_to_bigWig(char *outputfile, char *output_prefix, char *gtfile);
void make_wig(RefGenome *g, TYPE_WIGARRAY *array, char *outputfile, char *name, int binsize, int binnum, int chr, PWfile_Type wtype, int showzero);

#endif /* _OUTPUTFILE_H_ */
