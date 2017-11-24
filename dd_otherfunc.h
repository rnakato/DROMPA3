/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_OTHERFUNC_H_
#define _DD_OTHERFUNC_H_

#include "drompa_gv.h"
#include "dd_gv.h"

void dd_counttags(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample);
void dd_compare_intensity(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample);
void dd_multici(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample);
void dd_compare_genebody(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample);
void dd_travelling_ratio(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample);
#endif /* _DD_OTHERFUNC_H_ */
