/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_STROKE_H_
#define _DD_STROKE_H_

#include "drompa_gv.h"
#include "dd_gv.h"

void draw_region(DrParam *p, DDParam *d, SamplePair *sample, RefGenome *g, int s, int e, char *prefix, int region_no, int chr);
void merge_pdf(char *prefix);
int calc_pageheight(DrParam *p, DDParam *d);
void genome_overlook(DrParam *p, DDParam *d, RefGenome *g);

#endif /* _DD_STROKE_H_ */
