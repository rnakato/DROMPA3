/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_PROFILE_H_
#define _DD_PROFILE_H_

#include "drompa_gv.h"
#include "dd_gv.h"

void add_profile(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample, int chr);
void show_profile(DrParam *p, DDParam *d, SamplePair *sample);

#endif /* _DD_PROFILE_H_ */
