/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_INIT_H_
#define _DD_INIT_H_

#include "seq.h"
#include "drompa_gv.h"
#include "dd_gv.h"

void dd_argv_init(int argc, char *argv[], DrParam *p, DDParam *d, SamplePair **sample, RefGenome *g);

#endif /* _DD_INIT_H_ */
