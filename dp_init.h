/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DP_INIT_H_
#define _DP_INIT_H_

#include "drompa_gv.h"
#include "seq.h"

void dp_argv_init(int argc, char *argv[], DrParam *p, SamplePair **sample, RefGenome *g);

#endif /* _DP_INIT_H_ */
