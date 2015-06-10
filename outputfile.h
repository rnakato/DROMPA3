/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _OUTPUTFILE_H_
#define _OUTPUTFILE_H_

#include "seq.h"
#include "common.h"
#include "drompa_gv.h"

void output_bindata(char *dir, char *prefix, RefGenome *g, TYPE_WIGARRAY *array, char *gtfile, int binsize, int binnum, int chr, PWfile_Type wtype);
void make_binary(TYPE_WIGARRAY *array, char *outputfile, int binnum);
void make_bedGraph(RefGenome *g, TYPE_WIGARRAY *array, char *outputfile, char *name, int binsize, int binnum, int chr);
void convert_bedGraph_to_bigWig(char *outputfile, char *output_prefix, char *gtfile);
void make_wig(RefGenome *g, TYPE_WIGARRAY *array, char *outputfile, char *name, int binsize, int binnum, int chr, PWfile_Type wtype);

#endif /* _OUTPUTFILE_H_ */
