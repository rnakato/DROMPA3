/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DR_READFILE_H_
#define _DR_READFILE_H_

#include "seq.h"
#include "drompa_gv.h"

SamplePair *scan_samplestr(DrParam *p, char **, char **, int);
void dr_calc_global_param(DrParam *p, SamplePair *sample, RefGenome *g);
void dr_read_wigdata(DrParam *, SamplePair *sample, Samplefile *, RefGenome *g, int chr);
void dr_delete_wigdata(Samplefile *s);

#endif /* _DR_READFILE_H_ */
